function picplot(image_file_p, world_file_p)

image_file = imread(image_file_p);
world_file = worldfileread(world_file_p, 'planar', size(image_file));

%for xx = 1:2
%    [world_file.XWorldLimits(xx), world_file.YWorldLimits(xx)] = ll2utm(world_file.YWorldLimits(xx), world_file.XWorldLimits(xx));
%end
%figure;
mp = mapshow(image_file, world_file);
for ll = 1:2
    [mp.XData(ll), mp.YData(ll)] = ll2utm(mp.YData(ll), mp.XData(ll));
end

%%