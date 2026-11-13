# default settings
module Defaults
	# a nested module does not see the enclosing module's bindings, so the
	# names used below are imported explicitly
	const D = 128       # default bond dimension
	const rtol = 1e-12  # default relative truncation tolerance
	const tolgauge = 1e-14
	verbosity = 1
end
