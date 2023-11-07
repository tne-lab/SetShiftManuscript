function [t, y] = wienerrng2(a, ter, z, v)
  s = 1;
  delta_t = 0.001;
  max_t = 20;
  sq_st = sqrt(delta_t) * s;
  y = nan(1,max_t/delta_t);
  ind = 1;
  y(ind) = z * a;
  t = 0;
  while y(ind) <= a && y(ind) >= 0 && t <= max_t
      ind = ind + 1;
      y(ind) = y(ind-1) + v * delta_t + sq_st * randn;
      t = t + delta_t;
  end
  t = (t + ter) * sign(y(ind));
  y = rmmissing(y);
end