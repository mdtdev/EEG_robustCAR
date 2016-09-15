
function eg


  in_nm = '../../data/MG56_Seizure5.mat';
  load( in_nm );

  data                                            = sz.ECoG.Data';    % data
  [ n_channels n_samples ]                        = size( data );

  n_min_contributing_channels                         = 8;
  dt                                                  = 1;  % seconds
  n_threads                                           = 7;  % can have more than 1 thread per core.  I have two.
  [ d_rcar        ref_est_rcar    nn_ref_est_rcar  ]  = rCAR( data, dt, n_min_contributing_channels, n_threads );

end
