using Revise
using MultiPeriodDistOPF
using Profile
using ProfileView

# Import submodules with aliases
import .MultiPeriodDistOPF.computeOutputs as CO
import .MultiPeriodDistOPF.DDP as DDP
import .MultiPeriodDistOPF.DDPLinear as DDPLinear
import .MultiPeriodDistOPF.functionRetriever as FR
import .MultiPeriodDistOPF.helperFunctions as HF
import .MultiPeriodDistOPF.ModelBuilder as MB
import .MultiPeriodDistOPF.Playbook_of_MPOPF as Playbook
import .MultiPeriodDistOPF.parseOpenDSSFiles as Parser
import .MultiPeriodDistOPF.Plotter as Plotter
import .MultiPeriodDistOPF.Exporter as Exporter
import .MultiPeriodDistOPF.openDSSValidator as Validator

using JuMP
using Parameters

Revise.revise()

Revise.track(MultiPeriodDistOPF.computeOutputs)
Revise.track(MultiPeriodDistOPF.DDP)
Revise.track(MultiPeriodDistOPF.DDPLinear)
Revise.track(MultiPeriodDistOPF.functionRetriever)
Revise.track(MultiPeriodDistOPF.helperFunctions)
Revise.track(MultiPeriodDistOPF.ModelBuilder)
Revise.track(MultiPeriodDistOPF.Playbook_of_MPOPF)
Revise.track(MultiPeriodDistOPF.parseOpenDSSFiles)
Revise.track(MultiPeriodDistOPF.Plotter)
Revise.track(MultiPeriodDistOPF.Exporter)
Revise.track(MultiPeriodDistOPF.openDSSValidator)