# création de la solution
function create_solution(a,b,nb_inter)
    n_total = nb_inter + 1   # incluant a et b
    h = (b-a)/nb_inter
    x_temp = [a:h:b;]
    xf = x_temp[2:end-1]
    return n_total,h,xf
end

# création du point initial (perturbation du point final)
function create_initial_point_convex(a,b,nb_inter)
    n_total,h,xf = create_solution(a,b,nb_inter)
    legere_modif = x -> x^2
    x0 = legere_modif.(xf)
    return n_total,h,x0,xf
end 

# création du modèle 
function create_convex_JuMP_Model(;
                                 nb_inter::Integer=100,
                                 a=0,
                                 b=1,
                                 _string="undefined")
    m = Model()    
    (n_total,h,x0,xf) = create_initial_point_convex(a,b,nb_inter)# création du point initial en fonction du nombre d'intervalle
    n = length(x0)
    # @show n n_total x0
    create_initial_point_convex(a,b,nb_inter)
    @variable(m, x[1:n])
    
    @NLobjective(m, Min, sum( (x[j-1] - x[j+1])^2/(8*h) for j in 2:Integer(n-1)) + (x[2]-a)^2/(8*h) + (b-x[n-1])^2/(8*h) + (x[1]-a)^2/(4*h) + (b-x[n])^2/(4*h) )
    evaluator = JuMP.NLPEvaluator(m)
    MathOptInterface.initialize(evaluator, [:ExprGraph, :Hess])
    obj = MathOptInterface.objective_expr(evaluator)
    vec_var = JuMP.all_variables(m)
    @show length(vec_var) length(x0) n_total
    JuMP.set_start_value.(vec_var, x0)
    JuMP_nlp = MathOptNLPModel(m, name=_string*string(n))
    SPS_nlp = PartiallySeparableSolvers.PartionnedNLPModel(JuMP_nlp)
    return (m, JuMP_nlp, SPS_nlp)
end



#= FIN DE FONCTIONS=#

# Définition du modèle
nombre_intervals = 10
n,h,x0,xf = create_initial_point_convex(0,1,nombre_intervals)
(m,JuMP_nlp, SPS_nlp) = create_convex_JuMP_Model(;nb_inter=nombre_intervals)

# Lancement des solvers
sps_pbfgs, ges_pbfgs = PartiallySeparableSolvers.s_a_PBFGS(JuMP_nlp)
sps_psr1, ges_psr1 = PartiallySeparableSolvers.s_a_PSR1(JuMP_nlp)
sps_pbs, ges_pbs = PartiallySeparableSolvers.s_a_PBS(JuMP_nlp)

# Récupération des solutions
x_pbfgs = ges_pbfgs.solution
x_psr1 = ges_psr1.solution
x_pbs = ges_pbs.solution
# @test ges_pbfgs.solution ≈ xf atol=1e-6
# @test ges_pbs.solution ≈ xf atol=1e-6
# @test ges_psr1.solution ≈ xf atol=1e-6


# Récupération des métrices
H_pbfgs = zeros(n,n)
H_pbsr1 = zeros(n,n)
H_pbs = zeros(n,n)
H_pbfgs = PartiallySeparableNLPModel.construct_full_Hessian(sps_pbfgs.sps, sps_pbfgs.tpl_B[Int(sps_pbfgs.index)])
H_psr1 = PartiallySeparableNLPModel.construct_full_Hessian(sps_psr1.sps, sps_psr1.tpl_B[Int(sps_psr1.index)])
H_pbs = PartiallySeparableNLPModel.construct_full_Hessian(sps_pbs.sps, sps_pbs.tpl_B[Int(sps_pbs.index)])

# interval | itération | ni = 2
# 5 | 5-6
# 10 | 8
# 100 | 16-25
# 200 | 21-24

function create_convex_separable_JuMP_Model(;
		n::Integer=10,
		a=0,
		b=1,
		_string="undefined")
	m = Model()
	@variable(m, x[1:n])
	mod(n/2,2) != 0 && error("n n'est pas un multiple de 2")
	@NLobjective(m, Min, sum( (5*x[2*j-1] + x[2*j])^2 for j in 1:Integer(n/2)) )
	evaluator = JuMP.NLPEvaluator(m)
	MathOptInterface.initialize(evaluator, [:ExprGraph, :Hess])
	obj = MathOptInterface.objective_expr(evaluator)
	vec_var = JuMP.all_variables(m)	
	x0 = rand(n)
	JuMP.set_start_value.(vec_var, x0)
	JuMP_nlp = MathOptNLPModel(m, name=_string*string(n))
	SPS_nlp = PartiallySeparableSolvers.PartionnedNLPModel(JuMP_nlp)
	return (m, JuMP_nlp, SPS_nlp)
end

n = 100
(m,JuMP_nlp, SPS_nlp) = create_convex_separable_JuMP_Model(;n=n)

# Lancement des solvers
sps_pbfgs, ges_pbfgs = PartiallySeparableSolvers.s_a_PBFGS(JuMP_nlp)
sps_psr1, ges_psr1 = PartiallySeparableSolvers.s_a_PSR1(JuMP_nlp)
sps_pbs, ges_pbs = PartiallySeparableSolvers.s_a_PBS(JuMP_nlp)

# Récupération des solutions
x_pbfgs = ges_pbfgs.solution
x_psr1 = ges_psr1.solution
x_pbs = ges_pbs.solution
# @test ges_pbfgs.solution ≈ xf atol=1e-6
# @test ges_pbs.solution ≈ xf atol=1e-6
# @test ges_psr1.solution ≈ xf atol=1e-6

# Building the matrices
H_pbfgs = zeros(n,n)
H_pbsr1 = zeros(n,n)
H_pbs = zeros(n,n)
H_pbfgs = PartiallySeparableNLPModel.construct_full_Hessian(sps_pbfgs.sps, sps_pbfgs.tpl_B[Int(sps_pbfgs.index)])
H_psr1 = PartiallySeparableNLPModel.construct_full_Hessian(sps_psr1.sps, sps_psr1.tpl_B[Int(sps_psr1.index)])
H_pbs = PartiallySeparableNLPModel.construct_full_Hessian(sps_pbs.sps, sps_pbs.tpl_B[Int(sps_pbs.index)])

