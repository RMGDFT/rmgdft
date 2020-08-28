 
load input.xyz
cmd.bg_color('white')


cmd.viewport(1800, 1200)
set stick_color, black
preset.ball_and_stick(selection='all', mode=1)
#alter name Au, vdw=1.5
set sphere_scale, 0.8, name Au
set sphere_scale, 0.3, name C
set sphere_scale, 0.3, name S
set sphere_scale, 0.2, name H
color gold, name Au
color yellow, name S
