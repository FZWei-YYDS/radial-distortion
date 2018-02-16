function result = convergence_ditortion_center()
% convergence of the distortion center when applying the algorithm in
% "relation_rd_ru_sift_2" using a grid of points as initial values
%
pair_nr = 1;  %number of image pair (from 1 to 22)
cc1 = linspace(1800,2200,41);
cc2 = linspace(1300,1700,41);

k = 1;
for i=1:length(cc1)
    for j=1:length(cc2)
        cc0 = [cc1(i), cc2(j)];
        [ru,rd,cc,rmse,sum_areas] = relation_rd_ru_sift_2(pair_nr,cc0);
        result(k).cc = cc;
        result(k).rmse = rmse;
        result(k).sum_areas = sum_areas;
        fprintf('i:%d/%d, j:%d/%d\n',i,length(cc1),j,length(cc2));
        k = k + 1;
    end
    save result_convergence_cc1 result pair_nr
end

% plot the graph of the trajectories
% load('result_convergence_cc1');
rng(1);
figure(1); hold on; grid on;
for i=1:length(result)
    cx = rand(1,3);
    x = result(i).cc(:,1);
    y = result(i).cc(:,2);
    plot(x,y,'-','color',cx);
    plot(x(1),y(1),'x',x(end),y(end),'o','color',cx);
end