#TransientAnalysis class definition
struct TransientAnalysis <: Analysis#{{{
end #}}}

function UpdateParameters(analysis::TransientAnalysis,parameters::Parameters,md::model) #{{{

	AddParam(parameters, md.constants.yts, ConstantsYtsEnum)
	AddParam(parameters, md.timestepping.start_time*md.constants.yts, TimeEnum)
	AddParam(parameters, md.timestepping.final_time*md.constants.yts, TimesteppingFinalTimeEnum)
	AddParam(parameters, md.timestepping.time_step*md.constants.yts,  TimesteppingTimeStepEnum)
	AddParam(parameters, md.transient.isstressbalance, TransientIsstressbalanceEnum)
	AddParam(parameters, md.transient.ismasstransport, TransientIsmasstransportEnum)

end#}}}
function Core(analysis::TransientAnalysis,femmodel::FemModel)# {{{

	step      = FindParam(Int64, femmodel.parameters, StepEnum)
	time      = FindParam(Float64, femmodel.parameters, TimeEnum)
	finaltime = FindParam(Float64, femmodel.parameters, TimesteppingFinalTimeEnum)
	yts       = FindParam(Float64, femmodel.parameters, ConstantsYtsEnum)
	dt        = FindParam(Float64, femmodel.parameters, TimesteppingTimeStepEnum)

	isstressbalance = FindParam(Bool, femmodel.parameters, TransientIsstressbalanceEnum)
   ismasstransport = FindParam(Bool, femmodel.parameters, TransientIsmasstransportEnum)

   while(time < finaltime - (yts*eps(Float64))) #make sure we run up to finaltime.

		time+=dt
		AddParam(femmodel.parameters, time, TimeEnum)
		AddParam(femmodel.parameters, step, StepEnum)
		println("iteration ", step, "/", Int(ceil((finaltime-time)/dt))+step," time [yr]: ", time/yts, " (time step: ",  dt/yts, " [yr])")

      if(isstressbalance) Core(StressbalanceAnalysis(), femmodel) end
      if(ismasstransport) Core(MasstransportAnalysis(), femmodel) end
		MigrateGroundinglinex(femmodel)

		step+=1
   end

	println("=======================================")
	println("   Simulation completed successfully"   )
	println("=======================================")

end #}}}
