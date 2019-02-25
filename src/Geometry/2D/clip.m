function value = clip(x, min, max)
   if( min > max ); value = x; return; end
   if( x  < min ) ; value = min; return; end
   if( x  > max ) ; value = max; return; end
   value = x;
end