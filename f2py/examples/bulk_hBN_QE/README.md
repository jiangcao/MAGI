# bulk hBN 

## QE calculation 

`pw.scf.in` pw.x input file

```
 &control
    calculation='scf',
    restart_mode='from_scratch',
    pseudo_dir = '/home/jiacao/pseudopot/ONCVPSP/sg15/',
    outdir = '/usr/scratch2/tortin14/jiacao/tmp/',
    prefix='hbn',
    verbosity = 'high' ,
 /
 &system
    ibrav = 4,
    a = 2.49626563922,
    c = 6.00187404045,
    nat= 4,
    ntyp= 2,
    ecutwfc = 60.0,
!    nbnd = 32,
    vdw_corr = 'grimme-d2',
/
 &electrons
      diago_full_acc=.true. ,
      diagonalization='david' ,
      mixing_mode = 'plain' ,
      mixing_beta = 0.4,
      conv_thr =  1.0d-16,
/
ATOMIC_SPECIES
  B  0.0  B_ONCV_PBE_sr.upf
  N  0.0  N_ONCV_PBE_sr.upf
ATOMIC_POSITIONS (crystal)
 B             0.3333333000        0.6666667000        0.5000000000
 B             0.6666667000        0.3333333000        0.0000000000
 N             0.6666667000        0.3333333000        0.5000000000
 N             0.3333333000        0.6666667000        0.0000000000
K_POINTS (automatic)
 11  11  5   0 0 0
```

`wannier90` parameters in `EPW` input file

```
  dis_win_max = 30.0
  dis_win_min = -7.0
  dis_froz_max= 12.0
  dis_froz_min= -7.0
  proj(1)     = 'B:p,dxy'   
  proj(2)     = 'N:p,dxy'   
  wdata(1) = 'bands_plot = .true.'
  wdata(2) = 'begin kpoint_path'
  wdata(3) = 'K 0.3333333  0.3333333  0.00   G 0.0 0.00 0.00 '
  wdata(4) = 'G 0.0 0.00 0.00    Z 0.00 0.00 0.50'
  wdata(5) = 'Z 0.00 0.00 0.50   Q 0.3333333  0.3333333  0.50'
  wdata(6) = 'Q 0.3333333  0.3333333  0.50  L  0.5  0.0  0.50'
  wdata(7) = 'L 0.5  0.0  0.50   G 0.0 0.00 0.00 '
  wdata(8) = 'G 0.0 0.00 0.00    M 0.50 0.00 0.00'
  wdata(9) = 'M 0.50 0.00 0.00   K 0.3333333  0.3333333  0.00'
  wdata(10) = 'end kpoint_path'
  wdata(11) = 'bands_plot_format = gnuplot'
  wdata(12) = 'guiding_centres = .true.'
  wdata(13) = 'dis_num_iter      = 500'
  wdata(14) = 'num_print_cycles  = 10'
  wdata(15) = 'dis_mix_ratio     = 1.0'
  ```
