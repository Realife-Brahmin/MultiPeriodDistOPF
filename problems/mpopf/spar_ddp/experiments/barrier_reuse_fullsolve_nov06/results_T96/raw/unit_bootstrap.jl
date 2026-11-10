# Load actual production scan algebra without importing the full power model.
const UNIT_ROOT=get(ENV,"VALUE_REPLAY_NUMERICAL_ROOT",normpath(joinpath(@__DIR__,"../../..")))
include(joinpath(UNIT_ROOT,"spar_ddp/examples/power_system/pint_scan.jl"))
const UNIT_SOURCE=Meta.parseall(read(joinpath(UNIT_ROOT,"spar_ddp/examples/power_system/paired_pass.jl"),String))
function load_scan_units()
 methods_loaded=0
 for ex in UNIT_SOURCE.args
  ex isa Expr||continue
  if ex.head==:struct&&ex.args[2]==:PairWork
   Core.eval(Main,ex)
  elseif ex.head==:(=)&&ex.args[1] isa Expr&&ex.args[1].head==:call&&ex.args[1].args[1]==:PairWork
   Core.eval(Main,ex)
  elseif ex.head==:function&&ex.args[1] isa Expr&&ex.args[1].head==:call&&ex.args[1].args[1]==:_pp_combine!
   Core.eval(Main,ex);methods_loaded+=1
  end
 end
 @assert methods_loaded==2
end
load_scan_units()
