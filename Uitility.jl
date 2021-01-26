using JuMP
using Ipopt
using Plots
using PyPlot

struct Element
       name::String
       festWert::Float64
       gUnten::Float64
       gOben::Float64
end

struct Rohstoff
    name::String
    maxAvailable::Float64
    price::Float64
    elements::Vector{Element}
end

optimSol=Model[]
locOptimSol=Model[]

function elementContents(rohst::Rohstoff,element::Element)
    gehalt=Float64(0)
    for elem in rohst.elements
        if elem.name == element.name
            gehalt = elem.festWert
        end
    end
    return gehalt
end

contains(a::Element,b::Dict{String,Rohstoff}) = rohstoffeElement(a,b)
function rohstoffeElement(element, rohstoffe)
    alleRohstoffe=Dict{String,Rohstoff}()
    for roh in values(rohstoffe)
        for elem in roh.elements
            if elem.name == element.name
                push!(alleRohstoffe,roh)
            end
        end
    end
    return alleRohstoffe
end

#push a Rohstoff to an Dictionyry {String, ::Rohstoff}
Base.push!(a::Dict{String,Rohstoff}, b::Rohstoff) = dictPush(a,b)
Base.push!(a::Dict{String,Element},b::Element) = dictPush(a,b)
function dictPush(dict, value)
    push!(dict,value.name=>value)
end

#show Rohstoff as it should be
function Base.show(io::IO, b::Rohstoff)
    compact = get(io, :compact, false)
        show(io, b.name)
end

#Show an Vector of Rohstoffe as it should be
function Base.show(io::IO, b::Vector{Rohstoff})
    compact = get(io, :compact, false)
    if !compact
        for roh in b
            println("Rohstoff ", roh.name)
        end
    end
end

function Base.show(io::IO, element::Element)
    compact = get(io, :compact, false)
    if !compact
        println("Element ", element.name, " Min:", element.gUnten , " Max:", element.gOben , " actualValue: ", element.festWert)
    end
end

function simulate!(mengeSchmelze::Float64)

        m = Model(Ipopt.Optimizer)
        rohstoffe = Dict{String,Rohstoff}()
        elemente = Dict{String,Element}()

        mats = [Rohstoff("Kohle",1500,1.3,[Element("C",100,0,0)]),
                Rohstoff("Silizum",1500.0,2.5,[Element("Si",75.0,0.0,0.0)]),
                Rohstoff("FeCr",1500.0,2.2,[Element("Cr",65.0,0.0,0.0),
                                            Element("C",1.0,0.0,0.0)]),
                Rohstoff("FeCr0.2",1500.0,3.2,[Element("Cr",59.0,0.0,0.0),
                                            Element("C",0.02,0.0,0.0)]),
                Rohstoff("FeMn",1500.0,1.2,[Element("Mn",60.0,0.0,0.0)])
        ]
        Schmelze = Rohstoff("Schmelze",3000, 2.3, [
                                         Element("C",2.8,0.0,0.0),
                                         Element("Si",1.8,0.0,0.0),
                                         Element("Mn",0.8,0.0,0.0),
                                         Element("Cr",13.0,0.0,0.0)]
                           )

        elems = [Element("C",0,3.0,3.2),
                Element("Si",0,2.0, 2.2),
                Element("Mn",0,1.0, 1.5),
                Element("Cr",0,14.0, 15.0)]

        for mat in mats
            push!(rohstoffe,mat)
        end

        for elem in elems
            push!(elemente,elem)
        end

        @variable(m, 0 <= material[i=keys(rohstoffe)] <= rohstoffe[i].maxAvailable)

        @variable(m, elemente[i].gUnten <= element[i=keys(elemente)] <= elemente[i].gOben)
        @variable(m, MSchmelze == 3000)
        @variable(m, MGesamtSchmelze <= mengeSchmelze)

        @constraint(m, sum(material[i] for i in keys(rohstoffe)) + MSchmelze <= MGesamtSchmelze)

        tmpMats=Dict{String,Rohstoff}()
        for (key, elem) in elemente
            tmpMats = contains(elem,rohstoffe)
            @constraint(m, (sum((material[i] * elementContents(rohstoffe[i],elem))
                           for i in keys(tmpMats)) + (MSchmelze * elementContents(Schmelze,elem)))== MGesamtSchmelze * element[key])
        end

        @objective(m, Min, sum(material[i] * rohstoffe[i].price for i in keys(rohstoffe)) + MGesamtSchmelze * Schmelze.price)

        print(m)

        status = optimize!(m)
        println(termination_status(m))

        if termination_status(m) == JuMP.MOI.OPTIMAL
            push!(optimSol,m)
            println("optimale Lösung gefunden")
        elseif termination_status(m) == JuMP.MOI.LOCALLY_SOLVED
            push!(locOptimSol,m)
            println("Lokales Optimum gefunden")
        end


        println("Objective value: ", objective_value(m))
        for key in keys(rohstoffe)
            println("$key = ", round(value(material[key]),digits=1))
        end

        for key in keys(elemente)
            println("$key = ", round(value(element[key]),digits=2))
        end

        println("Menge Einlaufschmelze = ", value(MSchmelze))
        println("Menge Gesamt = ", value(MGesamtSchmelze))
        return m
end

function simulateMultipleMeltingWeigths(lowerBound::Int64,upperBound::Int64)
    objectives=Dict{Float64,Model}()
    for i in (lowerBound+10):upperBound
        ret = simulate!(Float64(i))
        push!(objectives, Float64(i)=>ret)
        if(termination_status(ret) == JuMP.MOI.OPTIMAL)
            break
        end
    end
    return objectives
end

function auswertungErgebnisse(results::Dict{Float64,Model})
    pyplot()
    x = []
    y = []

    minValue = 300000.00
    objValue=Float64
    minValueModel=Model()
    for (key, val) in results
        push!(x,key)
        push!(y,objective_value(val))
    end

    allVar = JuMP.all_variables(minValueModel)
    for i in allVar
        println(i, ":",value(i))
    end

    Plots.plot(x,y)
end

# results = simulateMultipleMeltingWeigths(3000,3550)
# auswertungErgebnisse(results)

simulate!(8000.0)
