using DataFrames
using CSV
using Statistics
using NLopt
using Zygote
using ArgParse

function load_data(data_file, mu_file)
    start_row = 3
    num_rows_per_read = 128
    dev_data = DataFrame!(CSV.File(data_file, header=1,
                    type=Float64, comment="CPS", ignoreemptylines=true, normalizenames=true));
    dev_data = dropmissing(dev_data)[:, 1:5];
    mu_data = DataFrame!(CSV.File(mu_file, header=1, datarow=3,
                type=Float64, ignoreemptylines=true, normalizenames=true));
    # mu_data = dropmissing(mu_data)[:, 1:2];
    return dev_data, mu_data
end

function compute_C1(μa, μs)
    return 3 * μa * μs
end

function compute_C2(μs, k0)
    return 6 * μs^2 * k0^2
end

function compute_K(τ, bfi, C1, C2)
    return √(C1+(C2*bfi*τ))
end

function compute_g2avg(df)
    col_names = [:Detector_1_Correlation, :Detector_2_Correlation, :Detector_3_Correlation, :Detector_4_Correlation]
    return map(mean, eachrow(df[:, col_names]))
end

function compute_g2fit(τ, β, bfi, r, C1, C2)
    K = compute_K(τ, bfi, C1, C2)
    return 1 + β * (exp(-r * K) / exp(-r * √C1))^2
end


function sqr_error(β, bfi, g2avg, τ, r, C1, C2)
    #bfi10 = 10^bfi
    bfi10 = bfi
    g2 = compute_g2fit(τ, β, bfi10, r, C1, C2)
    error = g2avg - g2
    return error^2
end

function sumsqrerror(x, grad, g2avg, τ, r, C1, C2)
    f(x) = sum(sqr_error.(x[1], x[2], g2avg, τ, r, C1, C2))
    grad .= (gradient(x) do x
            Zygote.forwarddiff(x) do x
                f(x)
            end
            end)[1]
    #grad .= gradient(f, x)[1]  # this was slow many more allocations
    return f(x)
end

function fitdata(tau, g2avg, μa, μs)
    k0 = 73919.8530152287
    r = 2.40
    
    C1 = compute_C1(μa, μs)
    C2 = compute_C2(μs, k0)  
    
    
    x = zeros(2)
    x .= [0.5, 1e-9]
    
    opt = Opt(:LD_MMA, length(x))
    
    f(x,g) = sumsqrerror(x, g, g2avg, tau, r, C1, C2)
    opt.lower_bounds = [0, 0]
    opt.upper_bounds = [1, 1e-6]
    
    opt.ftol_abs = 1e-15
    
    opt.min_objective = f
    (minf, minx, ret) = optimize(opt, x)
    numevals = opt.numevals  # number of function evaluations
    #println("Sum squared error $minf after $numevals iterations (returned $ret)")
    #println("beta: $(minx[1]) bfi: $(minx[2])")
    return minf, minx
end

function get_patient_params(dev_data, mu_data, num_rows_per_read)
    results = []
    for i in 1:convert(Int64, size(dev_data,1) / num_rows_per_read)
        μa, μs = mu_data[i, :mA], mu_data[i, :mS]
        start = 1 + (i-1) * num_rows_per_read
        last = i * num_rows_per_read

        g2avg = compute_g2avg(dev_data[start:last, :])
        tau = dev_data[start:last, :Time];
        minf, minx = fitdata(tau, g2avg, μa, μs)
        push!(results, [minf, minx[1], minx[2]])
    end
    results_mat = hcat(results...)'
    return results_mat
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--data", "-d"
            help = "file path containing g2 data"
            arg_type = String
            required = true
        "--mu", "-m"
            help = "file path containing muA and muS data"
            arg_type = String
            required = true
        "--out", "-o"
            help = "file to store results in"
            arg_type = String
            required = true
    end

    parsed_args = parse_args(ARGS, s)
    data_file = parsed_args["data"]
    mu_file = parsed_args["mu"]
    output_file = parsed_args["out"]
    
    num_rows_per_read = 128

    dev_data, mu_data = load_data(data_file, mu_file);

    results_mat = get_patient_params(dev_data, mu_data, num_rows_per_read);  # takes around 3 seconds to do 2446 measurements

    dfres = DataFrame(SSE = results_mat[:,1], Beta=results_mat[:,2], bfi=results_mat[:, 3])
    CSV.write(output_file, dfres)
end

main()



