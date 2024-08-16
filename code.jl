# Load necessary packages 
using NetSurvival, RateTables
using DataFrames
using Plots
using CategoricalArrays
using LaTeXStrings

# Age and follow-up distributions
age_fut = plot(
    histogram(colrec.age./365.241, label="Age"),
    histogram(colrec.time./365.241, label="Follow-up time"),
    size=(800,350)
)
savefig(age_fut, "colrec.png")

# Pohar Perme estimator
pp = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)

# Plot with confidence interval
conf_int = confint(pp; level = 0.05)
lower_bounds = [lower[1] for lower in conf_int]
upper_bounds = [upper[2] for upper in conf_int]
plot(pp.grid, pp.Sₑ, ribbon=(pp.Sₑ - lower_bounds, upper_bounds - pp.Sₑ), label = false, size=(800,400))
savefig("poharperme.png")

# nessie ess
breaks = [0, 45:5:90..., Inf]
colrec.agegr = cut(colrec.age./365.241, breaks)
ess, ess = nessie(@formula(Surv(time,status)~agegr), colrec, slopop)

# Crude mortality
crude_mortality = CrudeMortality(pp)

# Plot crude mortality
plot(pp.grid, crude_mortality.Mₑ, label = L"$\hat{M}_E$ Excess Mortality Rate", size=(800,400), ylims=(0.0,1.0))
plot!(pp.grid, crude_mortality.Mₚ, label = L"$\hat{M}_P$ Population Mortality Rate")
savefig("crude_mortality.png")

# New time and status 
colrec.time4000 .= 0.0
colrec.status4000 .= Bool(true)
for i in 1:nrow(colrec)
    colrec.time4000[i] = min(colrec.time[i], round(4000.0))
    if colrec.time[i] > 4000.0
        colrec.status4000[i] = false
    end
end

# New Pohar Perme
pp4000 = fit(PoharPerme, @formula(Surv(time4000,status4000)~1), colrec, slopop)
conf_int = confint(pp4000; level = 0.05)
lower_bounds = [lower[1] for lower in conf_int]
upper_bounds = [upper[2] for upper in conf_int]

crude_mortality4000 = CrudeMortality(pp4000)

p1 = plot(pp4000.grid, pp4000.Sₑ, ribbon=(pp4000.Sₑ - lower_bounds, upper_bounds - pp4000.Sₑ), label = false)
p2 = plot(pp4000.grid, crude_mortality4000.Mₑ, label = L"$\hat{M}_E$ Excess Mortality Rate")
p2 = plot!(pp4000.grid, crude_mortality4000.Mₚ, label = L"$\hat{M}_P$ Population Mortality Rate")
p = plot(p1,p2,size =(800,400), ylims=(0.0,1.0))
savefig(p, "new_plot.png")

# Grafféo test
fit(GraffeoTest, @formula(Surv(time4000,status4000)~sex), colrec, slopop)
fit(GraffeoTest, @formula(Surv(time4000,status4000)~site), colrec, slopop)
fit(GraffeoTest, @formula(Surv(time4000,status4000)~stage), colrec, slopop)
fit(GraffeoTest, @formula(Surv(time4000,status4000)~stage+Strata(sex)), colrec, slopop)

# Graphs comparing net survival
pp_sex = fit(PoharPerme, @formula(Surv(time4000,status4000)~sex), colrec, slopop)
male = pp_sex.estimator[1].Sₑ
conf_int_male = confint(pp_sex.estimator[1]; level = 0.05)
lower_bounds_male = [lower[1] for lower in conf_int_male]
upper_bounds_male = [upper[2] for upper in conf_int_male]

female =  pp_sex.estimator[2].Sₑ
conf_int_female = confint(pp_sex.estimator[2]; level = 0.05)
lower_bounds_female = [lower[1] for lower in conf_int_female]
upper_bounds_female = [upper[2] for upper in conf_int_female]

pp_sex_plot = plot(male, ribbon=(male - lower_bounds_male, upper_bounds_male - male), label="Male")
pp_sex_plot = plot!(female, ribbon=(female - lower_bounds_female, upper_bounds_female - female), label="Female")

pp_stage = fit(PoharPerme, @formula(Surv(time4000,status4000)~stage), colrec, slopop)
stage1 = pp_stage.estimator[1].Sₑ
conf_int1 = confint(pp_stage.estimator[1]; level = 0.05)
lower_bounds1 = [lower[1] for lower in conf_int1]
upper_bounds1 = [upper[2] for upper in conf_int1]

stage2 = pp_stage.estimator[2].Sₑ
conf_int2 = confint(pp_stage.estimator[2]; level = 0.05)
lower_bounds2 = [lower[1] for lower in conf_int2]
upper_bounds2 = [upper[2] for upper in conf_int2]

stage3 = pp_stage.estimator[3].Sₑ
conf_int3 = confint(pp_stage.estimator[3]; level = 0.05)
lower_bounds3 = [lower[1] for lower in conf_int3]
upper_bounds3 = [upper[2] for upper in conf_int3]

stage99 = pp_stage.estimator[4].Sₑ
conf_int99 = confint(pp_stage.estimator[4]; level = 0.05)
lower_bounds99 = [lower[1] for lower in conf_int99]
upper_bounds99 = [upper[2] for upper in conf_int99]

pp_stage_plot = plot(stage1, ribbon=(stage1 - lower_bounds1, upper_bounds1 - stage1), label="Stage 1")
pp_stage_plot = plot!(stage2, ribbon=(stage2 - lower_bounds2, upper_bounds2 - stage2), label="Stage 2")
pp_stage_plot = plot!(stage3, ribbon=(stage3 - lower_bounds3, upper_bounds3 - stage3), label="Stage 3")
pp_stage_plot = plot!(stage99, ribbon=(stage99 - lower_bounds99, upper_bounds99 - stage99), label="Stage 99")

plot(pp_sex_plot, pp_stage_plot, ylims=(0.0,1.0), size=(800,400))
savefig("net_comp.png")