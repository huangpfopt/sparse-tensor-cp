module batch_run

using CSV, DataFrames
using TSSOS
using Random
Random.seed!(1234);

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"matrix_IO.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"extract_atoms.jl")
include(proj_dir*"my_extract_atoms.jl")

tssos_dir = dirname(dirname(@__FILE__))*"\\term_sparsity\\"
include(tssos_dir*"tssos_model.jl")


export batch_clique_save,
       batch_comp_sparse_tssos_save 

function batch_clique_save(tens_dict, save_dir="results\\random\\clique\\")
    for name in sort([keys(tens_dict)...]) # name = "R_ex001_n6_d4_p0.4__"
        A = tens_dict[name]
        clique_save_path = save_dir*"checkClique_$(name).txt"
        n = parse(Int64,split(name,'_')[3][2:end])
        d = parse(Int64,split(name,'_')[4][2:end])
        _ = check_by_clique(A,n,d,clique_save_path)
    end
    return nothing
end
function check_by_clique(A,n,d,save_path::String)
    output, s = capture_solver_output(check_by_clique,(A,n,d))
    write_solver_output(s, save_path)
    return output
end
function check_by_clique(A,n,d)
    cpu = @elapsed begin
        mc = moments.get_maximal_cliques(A,n,d)
        might_CP = all(a->!isempty(a),[filter(c->issubset(item,c),mc) for item in moments.get_nonzero_entries(A)])
        end
    println("mc: ",mc)
    println("might_CP: ",might_CP)
    println("cpu: ", cpu)
    return mc, might_CP, cpu
end

function batch_comp_sparse_tssos_save(tens_dict, t, save_dir="results\\literature\\")
    for name in sort([keys(tens_dict)...])
        A = tens_dict[name]
        n = parse(Int64,split(name,'_')[2][2:end])
        d = parse(Int64,split(name,'_')[3][2:end])

        t0 = Int(ceil((d+1)/2))
        len = binomial(t0+n,t0)
        G = rand(len,len) 
        #G = randn(len,len)
        F = G'*G
        lvt = t0+t

        !isdir(save_dir*"t$(lvt)\\") ? mkdir(save_dir*"t$(lvt)\\") : 0

        F_dir  = save_dir*"t$(lvt)\\F\\"
        !isdir(F_dir) ? mkdir(F_dir) : 0
        F_save_path = F_dir*"F_$(name)_$(lvt).csv"
        matrix_IO.save_mat(F,F_save_path)
        
        moments_dir = save_dir*"t$(lvt)\\moments\\"
        !isdir(moments_dir) ? mkdir(moments_dir) : 0

        comp_sparse_save_path = save_dir*"t$(lvt)\\"*"ξₜᶜᵖ_$(name)_sp_$(lvt).txt"
        output = get_ξₜᶜᵖ(A,n,d,lvt,F,"sp",comp_sparse_save_path)
        mom_sparse_save_path = moments_dir*"Mom_$(name)_sp_$(lvt).csv"
        matrix_IO.save_moments(output[2],mom_sparse_save_path)

        comp_tssos_save_path = save_dir*"t$(lvt)\\"*"γₜᵗˢˢᵒˢ_$(name)_tssos_$(lvt).txt"
        output = get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F,comp_tssos_save_path) # get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F)
        mom_tssos_save_path = moments_dir*"Mom_$(name)_tssos_$(lvt).csv"
        #matrix_IO.save_mat(hcat(output[2]),mom_tssos_save_path) # Note that the output of tssos model is a moments vector
        matrix_IO.save_mat(TSSOS.get_moment_matrix(output[2], output[3])[1],mom_tssos_save_path) #Note that the output moment matrix has a different order of elements with the ones of ideal model, it need info1.basis to find the corresponding place, and TSSOS.get_nbasis is not Lexicographical ordered
    end
    return nothing
end
function get_ξₜᶜᵖ(A,n,d,t,F,flavour,save_path::String)
    output, s = capture_solver_output(get_ξₜᶜᵖ,(A,n,d,t,F,flavour)) # output = (ξₜᶜᵖ, mom_mat, cpu)
    write_solver_output(s, save_path)
    return output
end
function get_ξₜᶜᵖ(A,n,d,t,F,flavour)
    cpu = @elapsed begin ξₜᶜᵖ, mom_mat = cp_model.get_ξₜᶜᵖ(A,n,d,t,F,flavour) end 
    println("cpu: ", cpu)
    return ξₜᶜᵖ, mom_mat, cpu
end
function get_γₜᵗˢˢᵒˢ(A,n,d,t,F,save_path::String)
    output, s = capture_solver_output(get_γₜᵗˢˢᵒˢ,(A,n,d,t,F)) # output = (γₜᵗˢˢᵒˢ, mom_mat_tssos, info1, cpu)
    write_solver_output(s, save_path)
    return output
end
function get_γₜᵗˢˢᵒˢ(A,n,d,t,F)
    cpu = @elapsed begin γₜᵗˢˢᵒˢ, mom_mat_tssos, info1 = tssos_model.get_γₜᵗˢˢᵒˢ(A,n,d,t,F) end 
    println("cpu: ", cpu)
    return γₜᵗˢˢᵒˢ, mom_mat_tssos, info1, cpu
end

function batch_comp_dense_save(tens_dict, t, save_dir="results\\literature\\")
    for name in sort([keys(tens_dict)...])
        A = tens_dict[name]
        n = parse(Int64,split(name,'_')[2][2:end])
        d = parse(Int64,split(name,'_')[3][2:end])

        t0 = Int(ceil((d+1)/2))
        #len = binomial(t0+n,t0)
        #G = rand(len,len)
        #F = G'*G
        
        F_list = matrix_IO.load_mats(save_dir*"t$(t0)\\F\\")
        F = Matrix(F_list["F_$(name)_$(t0)"])
        lvt = t0+t

        #F_dir  = save_dir*"t$(lvt)\\F\\"
        #!isdir(F_dir) ? mkdir(F_dir) : 0
        #F_save_path = F_dir*"F_$(name)_$(lvt).csv"
        #matrix_IO.save_mat(F,F_save_path)
        !isdir(save_dir*"t$(lvt)\\") ? mkdir(save_dir*"t$(lvt)\\") : 0
        moments_dir = save_dir*"t$(lvt)\\moments\\"
        !isdir(moments_dir) ? mkdir(moments_dir) : 0

        comp_dense_save_path = save_dir*"t$(lvt)\\"*"ξₜᶜᵖ_$(name)_id_$(lvt).txt"
        output = get_ξₜᶜᵖ(A,n,d,lvt,F,"id",comp_dense_save_path)
        mom_dense_save_path = moments_dir*"Mom_$(name)_id_$(lvt).csv"
        matrix_IO.save_moments(output[2],mom_dense_save_path)
    end
    return nothing
end

function batch_extract_atoms(tens_dict,moments_dir) # "....\\moments\\"  
    fs = readdir(moments_dir)
    tfns = fs[map(f->(contains(f,"sp") && contains(f,".csv")),fs)]
    for tf in tfns
        spl = split(tf,'_')
        name = join(spl[2:end-2],"_") # example name
        n = parse(Int64,spl[3][2:end])
        d = parse(Int64,spl[4][2:end])
        model = spl[end-1]
        lvt = spl[end][1]

        A = tens_dict[name]
        # read file tf to load mom_mat
        mom_mat = matrix_IO.load_moments(moments_dir*tf)

        atoms_save_path = moments_dir*"Atom_$(name)_$(model)_$(lvt).txt"
        _ = process_moments(A,mom_mat,n,d,atoms_save_path)
    end
    return nothing
end

function process_moments(A,mom_mat,n,d,save_path::String)
    output, s = capture_solver_output(process_moments,(A,mom_mat,n,d)) # output = (mc, check, cpu)
    write_solver_output(s, save_path)
    return output
end
function process_moments(A,mom_mat,n,d)
    tol = 1e-6;
    mc = moments.get_maximal_cliques(A,n,d)
    println("mc: ",mc)
    cpu = @elapsed begin 
        n_sp, t_sp, mom_vals_sp = my_extract_atoms.proc_mom(mom_mat)
        check = my_extract_atoms.rank_check(n_sp, t_sp, d, mom_vals_sp, tol)
        try
            println("----my extract-----")
            ext_atoms_sp            = my_extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, tol)
            #ext_atoms_sp            = extract_atoms.get_atoms_robust(n_sp, t_sp, mom_vals_sp, tol, true)
            
            cents_sp, weights_sp    = my_extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
            A_ext_sp                = my_extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)
        
            error_sp = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
            println("recover_error: ", error_sp)
            println("cents_sp: ",cents_sp)
            println("weights_sp: ",weights_sp)
        catch
            println("No atom could be extracted")
        end
    end 

    cpu_hom = @elapsed begin 
        n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(mom_mat)
        try
            println("----hom extract-----")
            ext_atoms_sp            = extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, tol, true)
            
            cents_sp, weights_sp    = extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
            A_ext_sp                = extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)
        
            error_sp = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
            println("hom_recover_error: ", error_sp)
            println("hom_cents_sp: ",cents_sp)
            println("hom_weights_sp: ",weights_sp)
        catch
            println("No atom could be extracted")
        end
    end 
    println("cpu_hom: ", cpu_hom)

    println("check_flat: ",check)
    println("cpu: ", cpu)
    return mc, check, cpu
end

function batch_extract_id_atoms(tens_dict,moments_dir) # "....\\moments\\"  
    fs = readdir(moments_dir)
    tfns = fs[map(f->(contains(f,"id") && contains(f,".csv")),fs)]
    for tf in tfns
        spl = split(tf,'_')
        name = join(spl[2:end-2],"_") # example name
        n = parse(Int64,spl[3][2:end])
        d = parse(Int64,spl[4][2:end])
        model = spl[end-1]
        lvt = spl[end][1]

        A = tens_dict[name]
        # read file tf to load mom_mat
        mom_mat = matrix_IO.load_moments(moments_dir*tf)

        atoms_save_path = moments_dir*"Atom_$(name)_$(model)_$(lvt).txt"
        _ = process_id_moments(A,mom_mat,n,d,atoms_save_path)
    end
    return nothing
end

function process_id_moments(A,mom_mat,n,d,save_path::String)
    output, s = capture_solver_output(process_id_moments,(A,mom_mat,n,d)) # output = (check, cpu)
    write_solver_output(s, save_path)
    return output
end
function process_id_moments(A,mom_mat,n,d)
    tol = 1e-6;
    cpu = @elapsed begin 
        n_id, t_id, mom_vals_id = my_extract_atoms.proc_mom(mom_mat)
        check = my_extract_atoms.rank_check(n_id, t_id, d, mom_vals_id, tol)
        try
            println("----my extract-----")
            ext_atoms_id            = my_extract_atoms.get_atoms(n_id, t_id, mom_vals_id, tol)
            
            cents_id, weights_id    = my_extract_atoms.ext_centers_weights(ext_atoms_id)
            A_ext_id                = my_extract_atoms.recon_ten_id(cents_id, weights_id, d)
        
            error_id = sum([abs(A_ext_id[key]-A[key]) for key in keys(A)])
            println("recover_error: ", error_id)
            println("cents_id: ",cents_id)
            println("weights_id: ",weights_id)
        catch
            println("No atom could be extracted")
        end
    end 

    cpu_hom = @elapsed begin 
        n_id, t_id, mom_vals_id = extract_atoms.proc_mom(mom_mat)
        try
            println("----hom extract-----")
            ext_atoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, tol, true)
            
            cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)
            A_ext_id                = extract_atoms.recon_ten_id(cents_id, weights_id, d)
        
            error_id = sum([abs(A_ext_id[key]-A[key]) for key in keys(A)])
            println("hom_recover_error: ", error_id)
            println("hom_cents_id: ",cents_id)
            println("hom_weights_id: ",weights_id)
        catch
            println("No atom could be extracted")
        end
    end 
    println("cpu_hom: ", cpu_hom)

    println("check_flat: ",check)
    println("cpu: ", cpu)
    return check, cpu
end


function capture_solver_output(func,args)
    #warm-up
    output = func(args...)

    original_stdout = stdout;
    (rd,_) = redirect_stdout();
    output = func(args...) 
    s = []
    for rl in eachline(rd)
        push!(s,rl)
        if contains(rl,"cpu:")
            break #we need "break" to stop eachline(rd)
        end
    end
    redirect_stdout(original_stdout);
    return output, s
end
function write_solver_output(s,save_path)
    touch(save_path)
    open(save_path,"w") do io 
        for line in s
            write(io, line*"\n")
        end
    end;
end

### Summarize data-------------------------------------------------------------
function make_summary_clique(assets_dir)
    file_loc = assets_dir*"clique_summary.csv"
    #create_summary_file(file_loc)
    touch(file_loc)
    open(file_loc,"w") do io 
        write(io, "name,n,d,p,isCP,cpu,mc_sz \n")
    end

    fs = readdir(assets_dir)
    tfns = fs[map(f->contains(f,".txt"),fs)]
    for tf in tfns # tf: checkClique_R_ex366_n14_d6_p0.4__.txt
        spl = split(tf,'_')
        name = spl[3]
        n = parse(Int64,spl[4][2:end])
        d = parse(Int64,spl[5][2:end])
        p = parse(Float64,spl[6][2:end])

        # read file tf
        s = open(assets_dir*tf) do file 
            read(file, String) # readin the whole file as a long string
        end
        spl = split(s,'\n')
        isCP = split(spl[2])[2]
        cpu = parse(Float64,split(spl[3])[2]) #round(parse(Float64,split(spl[3])[2]), digits=4)

        mc = eval.(Meta.parse.(split(spl[1],':')[2]))
        mc_sz = [(sz, count(==(sz),length.(mc)) ) for sz in Set(length.(mc))] #[(1,2),(2,9)] means that has 2 cliques with 1 element and 9 cliques with 2 elements

        open(file_loc,"a") do io 
            write(io, "$name,$n,$d,$p,$isCP,$cpu,$mc_sz \n")
        end
    end
end

function make_summary_compare(assets_dir)
    file_loc = assets_dir*assets_dir[end-2:end-1]*"summary.csv"
    #create_summary_file(file_loc)
    touch(file_loc)
    open(file_loc,"w") do io 
        write(io, "name,n,d,model,obj,sdp_time,cpu,feas \n")
    end

    fs = readdir(assets_dir)
    tfns = fs[map(f->contains(f,".txt"),fs)]
    for tf in tfns #tf: ξₜᶜᵖ_ex01_n3_d3_17Fan-p7__id_2.txt
        spl = split(tf,'_')
        name = spl[2]
        n = parse(Int64,spl[3][2:end])
        d = parse(Int64,spl[4][2:end])
        model = spl[end-1]

        # read file tf
        s = open(assets_dir*tf) do file 
            read(file, String)
        end
        spl = split(s,'\n')[end-6:end-1]
        obj = round(parse(Float64,split(spl[end-1],':')[2]),digits=4)
        sdp_time = parse(Float64,split(spl[1],':')[2])
        cpu = parse(Float64,split(spl[end],':')[2]) #round(parse(Float64,split(spl[3])[2]), digits=4)
        feas = contains(spl[end-3],"Primal: FEASIBLE_POINT") && contains(spl[end-2],"Dual: FEASIBLE_POINT")
        
        open(file_loc,"a") do io 
            write(io, "$name,$n,$d,$model,$obj,$sdp_time,$cpu,$feas \n")
        end
    end
end



end