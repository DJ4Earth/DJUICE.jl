clear
close all
steps = [2];
clustername = 'totten';
saveflag = 1;
%Cluster parameters{{{
if strcmpi(clustername, 'andes')
   cluster=andes('numnodes',1,'cpuspernode',64, 'memory', 32);
   cluster.time = jobTime;
   waitonlock = 0;
elseif strcmpi(clustername, 'frontera')
   %cluster=frontera('numnodes', 1,'cpuspernode',56,'queue','flex');
   cluster=frontera('numnodes',3,'cpuspernode',56);
   cluster.time = jobTime;
   waitonlock = 0;
else
   cluster=generic('name',oshostname(),'np', 64);
   waitonlock = Inf;
end
clear clustername
org=organizer('repository',['./Models'],'prefix',['Amundsen_'],'steps',steps); clear steps;
%}}}
if perform(org, 'Sensitivity')% {{{
	md = loadmodel(org, 'Controls.mat');
	saveasstruct(md, './Models/Amundsen_Controls_dJUICE.mat');

	md.timestepping = timestepping();
	md.timestepping.start_time = 1995;
	md.timestepping.final_time = 1995.5;
	md.timestepping.time_step = 0.1;

	% manually set basal forcings
	md.basalforcings = linearbasalforcings();
	md.basalforcings.deepwater_melting_rate = 50.0;
	md.basalforcings.upperwater_melting_rate = 0.0;
	md.basalforcings.deepwater_elevation = -500.0;
	md.basalforcings.upperwater_elevation = 0.0;
	md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);
	md.basalforcings.perturbation_melting_rate = ones(md.mesh.numberofvertices, 1);
	% change the name of the first independent
	md.autodiff.independents{1}.name = 'BasalforcingsPerturbationMeltingRate';
	%md.groundingline.migration = 'None';

	md.toolkits=toolkits;
	md.cluster=generic('name',oshostname,'np',60);
	md.settings.output_frequency = 1;

	md=solve(md,'tr');

	savemodel(org,md);
end % }}}
if perform(org, 'Compare_with_DJUICE')% {{{
	md = loadmodel(org,['Sensitivity']);

	djuice = load('./Models/Amudsen_djuice_gradient.mat');

	grad_issm = rescalegradient(md ,md.results.TransientSolution(1).Gradient1);
	grad_djuice = rescalegradient(md, djuice.Gradient1);

	figure('Position',[0,400,1200,400])
   plotmodel(md, 'data', grad_issm,...
      'data', grad_djuice, ...
      'data', grad_issm - grad_djuice,...
      'colormap#3','bluewhitered',...
      'caxis#3',[-2.5e4,2.5e4],'caxis#1,2',[-1e5, 1e5],...
      'title', 'ISSM Sensitivity', 'title','DJUICE Sensitivity', 'title','ISSM - DJUICE',...
      'nlines', 1,...
      'mask#all', (md.mask.ocean_levelset<0 & md.mask.ice_levelset<0),...
      'axis#all', [-1.8275   -1.4059   -0.8326   -0.2671]*1e6,...
      'xlabel#all','x (m)','ylabel#all','y (m)')

   set(gcf,'color','w');

	if saveflag
		export_fig('ISSM_vs_DJUICE_melt.pdf')
	end
end % }}}
