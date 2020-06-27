using Random
using Statistics
using PyPlot
using SparseArrays

mutable struct Chain{T<:AbstractFloat}
    L::Int          # number of sites
    t::T            # in site terms of the hamiltonian (H_ii terms) 
    state::Int      # spatial state of the particle (position)
    V::T            # nn interaction terms from the hamiltonian (H_i,i+1 H_i,i-1 terms)
    H::Array{}      # hamiltonian
    g::Float64		#variational parameter
end

function δ(i::Int64,j::Int64)
    return i==j
end


function Hamiltonian(L::Int64,V::T,t::T) where D where T<:AbstractFloat #costruisce l'hamiltoniana in rappresentaz dei modi spaziali ; ok controllato
	#L=prod(I)
	H=zeros(L,L)
	V_=V/t
	for i in 1:L
	    for j in 1:L
	        H[i,j] = -(δ(i,j+1)+δ(i,j-1))+2V_*(j-1)*δ(i,j)
	    end
	end
	return H
end


function Chain(L::Int64,V::T,t::T,x0::Int64,g::Float64) where T<:AbstractFloat
	#costruttore dell'istanza Chain
	H=Hamiltonian(L,V,t)
	return Chain(L,t,x0,V,H,g) 
end

#=function H_l_x(x::Chain) 
    #local estimator of the energy
    g=x.g
    return sum(x.H[x.state,:].*exp.(-x.g * [i for i in 1:x.L])) / (exp(-g*x.state)) #ok controllato		
end =#

function H_l_x(x::Chain) 
    #local estimator of the energy
    g=x.g
    return sum(x.H[x.state,:].*exp.(-x.g * [i for i in 1:x.L])) / (exp(-g*x.state)) #ok controllato     
end 


function PdivP(x_::Int,x::Int,g::Float64) 
    #proposal ratio
	return (exp(-2*g*(x_-x)))
end

function T(x_::Int,x::Int,c::Chain) #ok controllato
    #trial 
    L=c.L
    x_==x+1 && x!=L && return 0.5
    x_==x-1 && x!=1 && return 0.5
    x_==x+1 && x==1 && return 1
    x_==x-1 && x==L && return 1
    return 0
end

function onemchop!(x::Chain) 
    # estraggo la config verso cui ho hopping con la trial definita sopra. il codice scritto è equivalente alla trial
	r=rand() 
	g=x.g
	state_=0 #il sito sul quale avrò l'hop, cioè la particella salta da x.state a state_ , per cui dev aggiornare x.state "salvandone" la cronologia SOLO se viene accettata la mossa
	if x.state!=1 && x.state!=L #non ai boundaries
		if r<0.5
			state_=x.state-1
            #println("Estratto +1")
		else 
			state_=x.state+1
            #println("Estratto -1")
		end
	else
		if x.state==1
			state_=2
            #printl("estratto +1 costretto")
		else 			#rimane solo il caso x.state==L
			state_=L-1
            #println("estratto -1 costretto")
		end
	end
	r=rand()
    A = min(1,PdivP(state_,x.state,g)*T(x.state,state_,x)/T(state_,x.state,x))
	if r<A 
		x.state=state_ #aggiorna lo stato compatibilmente con l'acceptance. Rispetta bilancio dettagliato
	end
end

mutable struct Measures{T<:AbstractFloat}
    it::Int                         # it contains the iteration step 
    ene::Vector{T}                  # a vector of length nmeas containing the energy measured at time it 
    pos::Vector{T}                  # a vector of length nmeas containing the position measured at time it 
    conf::Vector{Array{Int64,1}}    # a vector of length nmeas containing the spin configuration measured at time it 
end #inutile ridefinite un costruttore non standard, sono tutti parametri indipendenti


#function Chain(L::Int64,V::T,t::T,x0::Int64,g::Float64) where T<:AbstractFloat
function mcchain(L::Int64, t::T,                # linear size of the lattice + in site interaction terms
				V::T ,g::T ;                    #nn interaction hamiltonian term and variational parameter
                x0::Int=rand([i for i in 1:L]), # initial configuration 
                nterm::Int=1,                   # number of mchop for thermalization
                nmeas::Int=1000,                # number of measurments
                nhop::Int=1                     # number of mchop between measurements
                ) where {T<:Real}
#prendiamola prima misura del sistema dopo nterm passi. Le successive misure sono prese ogni nhop passi.In totale prendiamo nmeas misure
	output=Measures(0,zeros(Float64,nmeas),zeros(Float64,nmeas),Vector{Array{Int64,1}}(undef,nmeas)) #l'oggetto in cui salviamo le misure Measures(0,β,zeros(Float64,nmeas),zeros(Float64,nmeas),Vector{BitArray{1}}(undef,nmeas))
	# function Chain(I::NTuple{D,Int},V::Vector{T},t::T,x0::Int,g::T) where D where T<:AbstractFloat
	x=Chain(L,V,t,x0,g) #costruita bene, provato
	for j in 0:nterm #prendiamo la prima misura dopo il tempo di termalizzazione 
		onemchop!(x)
	end
	measure!(x,output) #effettua la prima misura
	for i in 2:nmeas #prendo nmeas misure, ciascuna a un tempo differente
		for j in 1:nhop #evolvo gli spin di nhop passi
			onemchop!(x)
		end
		measure!(x,output) #effettua la misura
	end
	return output , x 
end

function energy(x::Chain)
	return H_l_x(x::Chain)
end

function measure!(x::Chain,out::Measures)
    #effettua la misura, i.e. va a modificare gli attributi dell'"oggetto" out 
	out.it+=1                  #è stato inizializzato a 0. ogni volta che richiamo measure! aumento t di 1
	out.ene[out.it]=energy(x)
	out.pos[out.it]=x.state
	out.conf[out.it]=zeros(x.L) #inizializzo il vettore di config |x> = |0 , 0, ... 0>
	out.conf[out.it][x.state] = 1 #metto a 1 solo il ket associato ala posiz della particella, i.e. x.state
end

#=function block_field(output::Measures,n::Int,field::Symbol=:ene) where T 
	#divide output.field in blocchi (ex: energie) con lunghezze=potenze di n.  
    measures=copy(getfield(output,field))
    #Array di misure (measures) da dividere in blocchi di potenze n. PRENDI lenngth(measures) POTENZA DI n
    N=length(measures)
    function block_indices(n::Int,N::Int)
        #indici di start per dividere il singolo vettore di N misure in blocchi secondo potenza n (n,n^2,n^3,...)
        #N%n != 0 && return println("Non divisibili") inutile
        ind=Int[]
        function start_index(n,N)
            N<n && return
            start=start_index(n,N/n)
            push!(ind,Int(N))
            return
        end
        start_index(n,N)
        return ind
    end
    v_mean=[mean(measures[ind[i]+1:ind[i+1]]) for i in 1:Int(log(n,N))-1] #vettore di valori medi delle misure measures nei blocchi successivi
    variance_mean=std(v_mean,corrected=false)/sqrt(length(v_mean))
    return v_mean , variance_mean
end =#


function block_field(output::Measures,Lbin::Int,field::Symbol=:ene) where T #aggiungiamo field generico. Field può essere :energia, posizione, configurazione
	#restituisce un vettore e uno scalare:
	#								  1)vettore di valori medi di field in ogni bin 
	#								  3)deviazione standard della media (scalare)
    measures=copy(getfield(output,field))
    #Array di misure (measures) da dividere in blocchi di lunghezza Lbin. PRENDI lenngth(measures) MULTIPLO DI n
    N=Int(length(measures))
    Nbin=Int(N/Lbin)
    function block_indices_t(Lbin::Int,N::Int)
    #indici di stop per dividere il singolo vettore di N misure in blocchi di lunghezza Lbin
    N%Lbin != 0 && return println("Non divisibili")
    ind=Int[]
    function end_index_t(Lbin,N)
        N<Lbin && return
        end_ind=end_index_t(Lbin,N-Lbin)
        push!(ind,Int(N))
        return
    end
    end_index_t(Lbin,N)
        
    return ind
    end
    ind=block_indices_t(Lbin,N)
    blocks=[measures[i-Lbin+1:i] for i in ind]
    v_mean=[mean(blocks[i]) for i in 1:length(blocks)] #vettore di valori medi delle misure measures nei blocchi successivi
    #variance_block=[std(blocks[i],corrected=false) for i in 1:Nbin]
    variance_mean=std(v_mean)/sqrt(Nbin)
    return v_mean,variance_mean #blocks, v_mean , variance_block , variance_mean
end 


function L_bin_vect(n::Int,N::Int)
        #indici di start per L_bin secondo potenza n (n,n^2,n^3,...)
        N%n != 0 && return println("Non divisibili")
        ind=Int[]
        function stop_index(n,N)
            N<n && return
            start=stop_index(n,N/n)
            push!(ind,Int(N))
            return
        end
        stop_index(n,N)
        return ind
end