function new_unitFx=throwOutNanUnits(unitFx)

j=1;
for i=1:length(unitFx)
   if ~isnan(unitFx(i).no_led_FR) & unitFx(i).isGoodUnit==1
       new_unitFx(j)=unitFx(i);
       j=j+1;
   end    
end