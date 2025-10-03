clear
close all
steps = [2];
clustername = 'totten';
saveflag = 0;
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
org=organizer('repository',['./Models'],'prefix',['PIG_'],'steps',steps); clear steps;
%}}}

if perform(org, 'Sensitivity')% {{{

	md = loadmodel(org,['Control_drag']);

	saveasstruct(md, './Models/PIG_Control_drag_dJUICE.mat')

	md.friction.coefficient = 30*ones(size(md.friction.coefficient));

	md.inversion.iscontrol = 1;
	md.inversion.cost_functions=[101];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,1);
	md.inversion.cost_functions_coefficients(:,1)=1;
	md.inversion.maxsteps = 1;
	md.inversion.maxiter=1;
	md.inversion.incomplete_adjoint=0;

	md=solve(md,'Stressbalance');
	%plotmodel(md, 'data', abs(md.results.StressbalanceSolution.Gradient1), 'caxis', [0,5e-3],'colormap','jet')

	savemodel(org,md);
end %}}}
if perform(org, 'Compare_with_DJUICE')% {{{
	md = loadmodel(org,['Sensitivity']);

	djuice = load('./Models/PIG_djuice_gradient.mat');

	figure('Position',[0,400,1200,400])
	plotmodel(md, 'data', (md.results.StressbalanceSolution.Gradient1), ...
		'data',djuice.Gradient1, ...
		'data', md.results.StressbalanceSolution.Gradient1-djuice.Gradient1,...
		'colormap#3','bluewhitered',...
		'caxis#3',[-1e-4,1e-4],'caxis#1,2',[-1e-2,1e-2],...
		'title', 'ISSM Sensitivity', 'title','DJUICE Sensitivity', 'title','ISSM - DJUICE',...
		'nlines', 1,...
		'xlabel#all','x','ylabel#all','y')

	set(gcf,'color','w');
	if saveflag
		export_fig('ISSM_vs_DJUICE_Stressbalance.pdf')
	end
end % }}}

