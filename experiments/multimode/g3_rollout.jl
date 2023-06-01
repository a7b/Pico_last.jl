using Pico
using LinearAlgebra
const TRANSMON_LEVELS = 2
const CAVITY_LEVELS = 4

function cavity_state(level)
    state = zeros(CAVITY_LEVELS)
    state[level + 1] = 1.
    return state
end
TRANSMON_G = [1; zeros(TRANSMON_LEVELS - 1)]
TRANSMON_E = [zeros(1); 1; zeros(TRANSMON_LEVELS - 2)]
ψ1 = [kron(TRANSMON_G, cavity_state(0)), kron(TRANSMON_E, cavity_state(0))]

#ψf = [(kron(TRANSMON_G, cavity_state(0)) + kron(TRANSMON_G, cavity_state(4)))/sqrt(2), kron(TRANSMON_G, cavity_state(2))]

g2 = kron(TRANSMON_G, cavity_state(2))

prob = load_problem("/home/aditya/Pico.jl/data/multimode/fixed_time/no_guess/problems/g0_to_g1_T_120_dt_22.0_R_0.1_iter_3000ubound_0.1_00000.jld2")
#println(prob.objective_terms)
#prob = load_problem("data/twoqubit/CNOT_nelson_paper_00011.jld2")
controls_matrix = reduce(hcat, jth_order_controls(prob.trajectory, prob.system, 0; d2pi = false))

#infidelity = iso_infidelity(final_statei(prob.trajectory, system, i = 4), ket_to_iso(ψf[4]))

#display(final_statei(prob.trajectory, system, i = 3))
dts = [prob.trajectory.Δt for i = 1:prob.trajectory.T-1]
xf_1 = rollout(ket_to_iso(ψ1[1]), controls_matrix, dts, prob.system; integrator = Integrators.tenth_order_pade)[:, end]
#xf_2 = rollout(ket_to_iso(ψ1[2]), controls_matrix, dts, prob.system; integrator = Integrators.fourth_order_pade)[:, end]
println(fidelity(xf_1, ket_to_iso(g2)))
#println(fidelity(xf_2, ket_to_iso(ψf[2])))
