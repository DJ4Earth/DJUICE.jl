clear
md = loadmodel('./ModelAllControls.mat');

md.timestepping = timestepping();
md.timestepping.start_time = 1995;
md.timestepping.final_time = 1995.5;
md.timestepping.time_step = 0.1;

% manually set basal forcings
md.basalforcings = linearbasalforcings();
md.basalforcings.deepwater_melting_rate = 50.0
md.basalforcings.upperwater_melting_rate = 0.0
md.basalforcings.deepwater_elevation = -500.0
md.basalforcings.upperwater_elevation = 0.0
md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1)
md.basalforcings.perturbation_melting_rate = ones(md.mesh.numberofvertices, 1)
% change the name of the first independent
md.autodiff.independents{1}.name = 'BasalforcingsPerturbationMeltingRate';
%md.groundingline.migration = 'None';


md.toolkits=toolkits;
md.cluster=generic('name',oshostname,'np',60);
md.settings.output_frequency = 1;

%md.inversion.iscontrol = 0;
%md.autodiff.isautodiff = 0;
%md.verbose = verbose('All');

md=solve(md,'tr');
save ./Model_AD.mat md
saveasstruct(md, 'issm_ad.mat');
