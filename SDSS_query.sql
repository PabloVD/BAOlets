SELECT ra,dec,z
FROM SpecPhoto

WHERE

class=’galaxy’
AND z>0.15 and z<0.43
AND ra<225 and ra>125
AND dec>0 and dec<60
AND ((cmodelmag_r < 13.5 + (0.7*(modelMag_g-modelMag_r) +
1.2*(( modelMag_r-modelMag_i)-0.18))/0.3
AND (modelMag_r-modelMag_i -(modelMag_g-modelMag_r)/4.0 - 0.18) >
-0.2 and (modelMag_r-modelMag_i -(modelMag_g-modelMag_r)/4.0 - 0.18) < 0.2
AND cmodelmag_r > 16.0 and cmodelmag_r < 19.6
AND psfmag_r-cmodelmag_r > 0.3
AND TILE >=10324.0
AND ZWARNING = 0)

OR

(cmodelmag_r < 13.5 + (0.7*(modelMag_g-modelMag_r) +
1.2*((modelMag_r-modelMag_i)-0.18))/0.3
AND (modelMag_r-modelMag_i -(modelMag_g-modelMag_r)/4.0 - 0.18) >
-0.2 and (modelMag_r-modelMag_i -(modelMag_g-modelMag_r)/4.0 - 0.18) < 0.2
AND cmodelmag_r > 16.0 and cmodelmag_r < 19.6
AND psfmag_r-cmodelmag_r > 0.3
AND Boss_target1=0
AND ZWARNING = 0))

ORDER BY z ASC
