using dJUICE
using MAT
using Enzyme

# params = dJUICE.Parameters(Dict{dJUICE.IssmEnum, dJUICE.Parameter}())
params = dJUICE.Parameters(
    Dict{dJUICE.IssmEnum, dJUICE.Parameter}(),
    Dict{dJUICE.IssmEnum,dJUICE.DoubleParam}()
)
dJUICE.AddParam(params, 10.0, dJUICE.BasalforcingsOceanSalinityEnum)

dparams = dJUICE.Parameters(
    Dict{dJUICE.IssmEnum, dJUICE.Parameter}(),
    Dict{dJUICE.IssmEnum,dJUICE.DoubleParam}()
)
dJUICE.AddParam(dparams, 0.0, dJUICE.BasalforcingsOceanSalinityEnum)


function f(params, x)
    return x * dJUICE.FindParam(Float64, params, dJUICE.BasalforcingsOceanSalinityEnum)
end

@show f(params, 5.0)

using Enzyme
@show autodiff(f, Active, Duplicated(params, dparams), Active(5.0))

@show dparams