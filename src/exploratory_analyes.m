clear 
close all
addpath('utilities')
project = 'eveBAC_2019-05-14';
load(['../dat/' project '/nucleus_struct.mat'])
figPath = '../fig/exploratory/';
mkdir(figPath)

%%% Add some useful flags
Tres = nucleus_struct(1).TresInterp;
for i = 1:numel(nucleus_struct)
    time = nucleus_struct(i).time;
    fluo = nucleus_struct(i).fluo;
            
    nucleus_struct(i).nc14_flag = max(time) - min(time)>= 20*60;
    nucleus_struct(i).n_on = sum(max(~isnan(fluo),[],2));
    start_vec = NaN(2,1);
    stop_vec = NaN(2,1);
    gap_vec = NaN(2,1);
    gap5_vec = NaN(2,1);
    for j = 1:2
        fv = fluo(j,:);
        if any(~isnan(fv)) 
            start_i = find(~isnan(fv),1);
            stop_i = find(~isnan(fv),1,'last');
            start_vec(j) = time(start_i);
            stop_vec(j) = time(stop_i);
            dur_filt = conv(~isnan(fv),ones(1,5),'same');
            gap =  max(diff(find(dur_filt>0)))*Tres;
            gap5 =  max(diff(find(dur_filt==5)))*Tres;
            if ~isempty(gap) && ~isempty(gap5)
                gap_vec(j) = gap;
                gap5_vec(j) = gap5;
            else
                gap_vec(j) = 0;
                gap5_vec(j) = 0;
            end
        end
    end
    nucleus_struct(i).start_time = start_vec;
    nucleus_struct(i).late_start = start_vec > 30*60;
    nucleus_struct(i).stop_time = stop_vec;
    nucleus_struct(i).early_stop = stop_vec <= 30*60;
    nucleus_struct(i).biggest_gap = gap_vec;
    nucleus_struct(i).big_gap = gap_vec >= 15*60&gap5_vec >= 15*60;
    
    % calculate nucleus AP position
    nucleus_struct(i).apMean = nanmean(nucleus_struct(i).apPos);
end
    

% get average fraction of nuclei with 0,1, and 2 active spots as a function
% of AP
ap_vec = 1:100;
active_ap_array = zeros(3,100);
nc_ap_vec = [nucleus_struct.apMean];
nc14_flag = [nucleus_struct.nc14_flag];
for a = 1:numel(ap_vec)
    ap_filter = round(nc_ap_vec*100)==ap_vec(a)&nc14_flag;
    active_ap_array(1,a) = sum([nucleus_struct(ap_filter).n_on]==0)/sum(ap_filter);
    active_ap_array(2,a) = sum([nucleus_struct(ap_filter).n_on]==1)/sum(ap_filter);
    active_ap_array(3,a) = sum([nucleus_struct(ap_filter).n_on]==2)/sum(ap_filter);
end

% number of active spots as a function of AP
n_on_fig = figure;
n_on_fig.Position = [100 100 1028 512];
cm = brewermap(8,'Set1');
colormap(cm)
hold on
% for i = 1:3
%     plot(ap_vec,active_ap_array(i,:),'LineWidth',1.3,'Color',cm(3-i+1,:))
% end
a = area(repmat(ap_vec',1,3),active_ap_array');
for i = 1:numel(a)
    a(i).FaceColor = cm(3-i+1,:);
    a(i).FaceAlpha = .7;
end
grid on
legend('0 on','1 on', '2 on','Location','northeast')
grid on
xlabel('% embryo length')
ylabel('share') 
xlim([15,60])
saveas(n_on_fig,[figPath 'n_on_fig.png'])
%%         
% calculate max times for each set
set_vec = [nucleus_struct.setID];
set_index = unique(set_vec);
max_time_vec = NaN(1,numel(set_index));
for i = 1:numel(set_index)
    max_time_vec(i) = nanmax(nanmax([nucleus_struct(set_vec==set_index(i)).time_interp]));
end
% correlation in on times amongst sister
two_spot_flag  =nc14_flag&[nucleus_struct.n_on]==2 & sum([nucleus_struct.inference_flag])==2;
ref_flag = sum([nucleus_struct.inference_flag])>0&nc14_flag;
two_spot_indices = find(two_spot_flag);
nc_x_vec = [nucleus_struct.xMean];
nc_y_vec = [nucleus_struct.yMean];
nc_ap_vec = [nucleus_struct.apMean];
setID_vec = [nucleus_struct.setID];
% initialize structures
sister_struct = struct;
nn_struct = struct;
rng(123);
bkg_iter = 1;
for i = two_spot_indices
    % sisters first
    setID = nucleus_struct(i).setID;
    max_time = max_time_vec(set_index==setID);
    fluo_array = nucleus_struct(i).fluo_interp;
    time_vec = nucleus_struct(i).time_interp(1,:);
    sister_on_vec = [];
    sister_off_vec = [];
    indices = randsample(1:2,2,false);
    for j = indices
        fluo = fluo_array(j,:);
        sister_on_vec = [sister_on_vec time_vec(find(~isnan(fluo),1))];
        off_time = time_vec(find(~isnan(fluo),1,'last'));      
        sister_off_vec = [sister_off_vec off_time];
    end
    sister_off_vec(sister_off_vec>=max_time) = NaN;
    sister_struct(bkg_iter).on_times = sister_on_vec;
    sister_struct(bkg_iter).off_times = sister_off_vec;
    sister_mf = nanmean(fluo_array,2)';
    sister_struct(bkg_iter).mean_fluo = sister_mf;
    sister_var = nanvar(fluo_array,[],2)';
    sister_struct(bkg_iter).var_fluo = sister_var;
    % compare to nearest neighbor
    sister_id = ceil(rand()*2);
    setID = setID_vec(i);
    x_vec = nc_x_vec - nc_x_vec(i);
    y_vec = nc_y_vec - nc_y_vec(i);
    r_vec = sqrt(x_vec.^2 + y_vec.^2);
    r_vec(i) = Inf;
    r_vec(setID_vec~=setID|~ref_flag) = Inf;
    % find closest nucleus
    [~, mi_gap] = min(r_vec);
    options = find(nucleus_struct(mi_gap).inference_flag==1);
    nn_ind = randsample(options,1);
    nn_time_vec = nucleus_struct(mi_gap).time_interp(1,:);
    nn_fluo = nucleus_struct(mi_gap).fluo_interp(nn_ind,:);
    % calculate matrics
    nn_on_vec = [sister_on_vec(sister_id) nn_time_vec(find(~isnan(nn_fluo),1))];
    nn_off_vec = [sister_off_vec(sister_id) nn_time_vec(find(~isnan(nn_fluo),1,'last'))];
    nn_off_vec(nn_off_vec>=max_time) = NaN;
    nn_struct(bkg_iter).on_times = nn_on_vec;
    nn_struct(bkg_iter).off_times = nn_off_vec;
    nn_struct(bkg_iter).mean_fluo = [sister_mf(sister_id) nanmean(nn_fluo)];
    nn_struct(bkg_iter).var_fluo = [sister_var(sister_id) nanvar(nn_fluo)];    
    bkg_iter = bkg_iter + 1;    
end   

%%% on and off times
t_axis = linspace(0,50);
% concatenate results
sister_on_array = vertcat(sister_struct.on_times);
sister_off_array = vertcat(sister_struct.off_times);
nn_on_array = vertcat(nn_struct.on_times);
nn_off_array = vertcat(nn_struct.off_times);
% linear fit
nn_on_vec1 = nn_on_array(:,1)/60;
nn_on_vec2 = nn_on_array(:,2)/60;
p = polyfit(nn_on_vec1,nn_on_vec2,1);
nn2_pd = polyval(p,t_axis);

nn_on_fig = figure;
hold on
scatter(nn_on_array(:,1)/60,nn_on_array(:,2)/60,30,'o','MarkerFaceColor',cm(2,:),'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0)
plot(t_axis,nn2_pd,'Color','black')
legend('nearest neighbor spots','linear fit','Location','northwest')
xlabel('on time (spot 1)')
ylabel('on time (spot 2)')
grid on
box on
set(gca,'Fontsize',12)
saveas(nn_on_fig,[figPath 'nn_on_time_scatter.png'])

% linear fit
s_on_vec1 = sister_on_array(:,1)/60;
s_on_vec2 = sister_on_array(:,2)/60;
p = polyfit(s_on_vec1,s_on_vec2,1);
s2_pd = polyval(p,t_axis);

sister_on_fig = figure;
hold on
scatter(s_on_vec1,s_on_vec2,30,'o','MarkerFaceColor',cm(1,:),'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0)
plot(t_axis,s2_pd,'Color','black')
legend('sister spots','linear fit','Location','northwest')
xlabel('on time (spot 1)')
ylabel('on time (spot 2)')
grid on
box on
set(gca,'Fontsize',12)
saveas(sister_on_fig,[figPath 'sister_on_time_scatter.png'])

% linear fit
nn_off_vec1 = nn_off_array(:,1)/60;
nn_off_vec2 = nn_off_array(:,2)/60;
nan_ft = ~isnan(nn_off_vec1)&~isnan(nn_off_vec2);
p = polyfit(nn_off_vec1(nan_ft),nn_off_vec2(nan_ft),1);
nn2_pd = polyval(p,t_axis);

nn_off_fig = figure;
hold on
scatter(nn_off_array(:,1)/60,nn_off_array(:,2)/60,30,'o','MarkerFaceColor',cm(2,:),'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0)
plot(t_axis,nn2_pd,'Color','black')
legend('nearest neighbor spots','linear fit','Location','northwest')
xlabel('off time (spot 1)')
ylabel('off time (spot 2)')
grid on
box on
set(gca,'Fontsize',12)
saveas(nn_off_fig,[figPath 'nn_off_time_scatter.png'])

% linear fit
s_off_vec1 = sister_off_array(:,1)/60;
s_off_vec2 = sister_off_array(:,2)/60;
nan_ft = ~isnan(s_off_vec1)&~isnan(s_off_vec2);
p = polyfit(s_off_vec1(nan_ft),s_off_vec2(nan_ft),1);
s2_pd = polyval(p,t_axis);

sister_off_fig = figure;
hold on
scatter(s_off_vec1,s_off_vec2,30,'o','MarkerFaceColor',cm(1,:),'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0)
plot(t_axis,s2_pd,'Color','black')
legend('sister spots','linear fit','Location','northwest')
xlabel('off time (spot 1)')
ylabel('off time (spot 2)')
grid on
box on
set(gca,'Fontsize',12)
saveas(sister_off_fig,[figPath 'sister_off_time_scatter.png'])


%%
particle_id_vec = [nucleus_struct.ParticleID];
particle_flag_vec = max(~isnan(particle_id_vec))>0;
gap_flag_vec = [nucleus_struct.big_gap];
[gap_flag_vec, mi_gap] = max(gap_flag_vec);
late_flag_vec = [nucleus_struct.late_start];
[late_flag_vec, mi_late] = max(late_flag_vec);

bkg_filter = ~late_flag_vec&~gap_flag_vec & particle_flag_vec & nc14_flag;
gap_filter = ~late_flag_vec & gap_flag_vec & particle_flag_vec & nc14_flag;
late_filter = late_flag_vec & ~gap_flag_vec & particle_flag_vec & nc14_flag;


% initialize and populate arrays to track single particle trajectories
bkg_ap_array = NaN(200,sum(bkg_filter));
bkg_time_array = NaN(200,sum(bkg_filter));
bkg_fluo_array = NaN(200,sum(bkg_filter));
gap_ap_array = NaN(200,sum(gap_filter));
gap_time_array = NaN(200,sum(gap_filter));
gap_fluo_array = NaN(200,sum(gap_filter));
late_ap_array = NaN(200,sum(gap_filter));
late_time_array = NaN(200,sum(gap_filter));
late_fluo_array = NaN(200,sum(gap_filter));
iter = 1;
for i = find(bkg_filter)
    fluo = nucleus_struct(i).fluo(mi_gap(i),:);
    time = nucleus_struct(i).time(mi_gap(i),:);
    ap = nucleus_struct(i).apPos(mi_gap(i),:);
    bkg_ap_array(1:numel(ap),iter) = ap;
    bkg_time_array(1:numel(ap),iter) = time;
    bkg_fluo_array(1:numel(ap),iter) = fluo;
    iter = iter + 1;
end
iter = 1;
for i = find(gap_filter)
    fluo = nucleus_struct(i).fluo(mi_gap(i),:);
    time = nucleus_struct(i).time;%(mi(i),:);
    ap = nucleus_struct(i).apPos;%(mi(i),:);
    gap_ap_array(1:numel(ap),iter) = ap;
    gap_time_array(1:numel(ap),iter) = time;
    gap_fluo_array(1:numel(ap),iter) = fluo;
    iter = iter + 1;
end
iter = 1;
for i = find(late_filter)
    fluo = nucleus_struct(i).fluo(mi_late(i),:);
    time = nucleus_struct(i).time;%(mi(i),:);
    ap = nucleus_struct(i).apPos;%(mi(i),:);
    late_ap_array(1:numel(ap),iter) = ap;
    late_time_array(1:numel(ap),iter) = time;
    late_fluo_array(1:numel(ap),iter) = fluo;
    iter = iter + 1;
end


close all
cm1 = brewermap(11,'Spectral');

base_figure = figure;
base_figure.Position = [100 100 1028 512];
hold on
% background 
p1 = plot(bkg_ap_array,bkg_time_array/60,'Color',[.5 .5 .5 .1],'LineWidth',1.5);
bkg_ap_array_spot = bkg_ap_array;
bkg_ap_array_spot(isnan(bkg_fluo_array)) = NaN;
bkg_time_array_spot = bkg_time_array;
bkg_time_array_spot(isnan(bkg_fluo_array)) = NaN;
p2 = plot(bkg_ap_array_spot,bkg_time_array_spot/60,'Color',[cm1(end-1,:) .3],'LineWidth',1.5);
set(gca, 'YDir','reverse')
xlim([.25 .62]);
ylim([5 50])
grid on
xlabel('% embryo length')
ylabel('minutes into nc14')
legend([p1(1) p2(1)],'nucleus tracks','active periods','Location','northwest')
set(gca,'Fontsize',12)
saveas(base_figure,[figPath 'base_trace_map.png'])

gap_figure = figure;
gap_figure.Position = [100 100 1028 512];
hold on
% background 
scatter(bkg_ap_array(:),bkg_time_array(:)/60,10,'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);%,'LineWidth',1.5)
bkg_ap_array_spot = bkg_ap_array;
bkg_ap_array_spot(isnan(bkg_fluo_array)) = NaN;
bkg_time_array_spot = bkg_time_array;
bkg_time_array_spot(isnan(bkg_fluo_array)) = NaN;
scatter(bkg_ap_array_spot(:),bkg_time_array_spot(:)/60,'MarkerFaceColor',[cm1(end-1,:)],'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);%,'LineWidth',1.5)
% gap traces
p1 = plot(gap_ap_array,gap_time_array/60,'Color',[cm1(end,:) .4],'LineWidth',2);
gap_ap_array_spot = gap_ap_array;
gap_ap_array_spot(isnan(gap_fluo_array)) = NaN;
gap_time_array_spot = gap_time_array;
gap_time_array_spot(isnan(gap_fluo_array)) = NaN;
p2 = plot(gap_ap_array_spot,gap_time_array_spot/60,'Color',[cm1(4,:) 1],'LineWidth',2);
% set axis characteristics
set(gca, 'YDir','reverse')
xlim([.25 .62]);
ylim([5 50])
grid on
xlabel('% embryo length')
ylabel('minutes into nc14')
legend([p1(1) p2(1)],'nucleus tracks','active periods (gap nuclei)','Location','northwest')
set(gca,'Fontsize',12)
saveas(gap_figure,[figPath 'gap_trace_map.png'])

late_figure = figure;
late_figure.Position = [100 100 1028 512];
hold on
% background 
scatter(bkg_ap_array(:),bkg_time_array(:)/60,10,'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);%,'LineWidth',1.5)
bkg_ap_array_spot = bkg_ap_array;
bkg_ap_array_spot(isnan(bkg_fluo_array)) = NaN;
bkg_time_array_spot = bkg_time_array;
bkg_time_array_spot(isnan(bkg_fluo_array)) = NaN;
scatter(bkg_ap_array_spot(:),bkg_time_array_spot(:)/60,'MarkerFaceColor',[cm1(end-1,:)],'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0);%,'LineWidth',1.5)
% gap traces
p1 = plot(late_ap_array,late_time_array/60,'Color',[cm1(end,:) .4],'LineWidth',2);
late_ap_array_spot = late_ap_array;
late_ap_array_spot(isnan(late_fluo_array)) = NaN;
late_time_array_spot = late_time_array;
late_time_array_spot(isnan(late_fluo_array)) = NaN;
p2 = plot(late_ap_array_spot,late_time_array_spot/60,'Color',[cm1(1,:) 1],'LineWidth',2);
% set axis characteristics
set(gca, 'YDir','reverse')
xlim([.25 .62]);
ylim([5 50])
grid on
xlabel('% embryo length')
ylabel('minutes into nc14')
legend([p1(1) p2(1)],'nucleus tracks','active periods (late nuclei)','Location','northwest')
set(gca,'Fontsize',12)
saveas(late_figure,[figPath 'late_trace_map.png'])

%% mean fluorescence
sister_mf_array = vertcat(sister_struct.mean_fluo);
nn_mf_array = vertcat(nn_struct.mean_fluo);

mean_fluo_fig = figure;
hold on
scatter(sister_mf_array(:,1),sister_mf_array(:,2),30,'o','MarkerFaceColor',cm(1,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0)
scatter(nn_mf_array(:,1),nn_mf_array(:,2),30,'o','MarkerFaceColor',cm(2,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0)
legend('sister spots','nearest neighbor spots','Location','northwest')
xlabel('mean activity (spot 1)')
ylabel('mean activity (spot 2)')
grid on
set(gca,'Fontsize',12)
saveas(mean_fluo_fig,[figPath 'sister_mean_fluo_scatter.png'])

sister_var_array = sqrt(vertcat(sister_struct.var_fluo));
nn_var_array = sqrt(vertcat(nn_struct.var_fluo));

var_fluo_fig = figure;
hold on
scatter(sister_var_array(:,1),sister_var_array(:,2),30,'o','MarkerFaceColor',cm(1,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0)
scatter(nn_var_array(:,1),nn_var_array(:,2),30,'o','MarkerFaceColor',cm(2,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0)
legend('sister spots','nearest neighbor spots','Location','northwest')
xlabel('variance (spot 1)')
ylabel('variance (spot 2)')
grid on
set(gca,'Fontsize',12)
saveas(mean_fluo_fig,[figPath 'sister_mean_fluo_scatter.png'])



