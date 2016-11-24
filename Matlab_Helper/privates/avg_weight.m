function weight = avg_weight(fvalue, vert)
avg_value = (fvalue(1:end-1)+fvalue(2:end))/2.0;
diff = vert(1:end-1, 1:3) - vert(2:end, 1:3);
dist = sqrt(sum(diff.*diff, 2));
sum_dist = sum(dist);

weight = sum(avg_value.*dist) / sum_dist;