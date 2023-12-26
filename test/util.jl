
f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D) = SpectrumFunction(ϵ->f(D, ϵ), lb=-D, ub=D)
spectrum_func() = spectrum_func(10)
