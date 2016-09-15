function [data_out, ref_est, nn_ref_est] = robustCAR4_winsz(data_in, dt, n_min_contributing_channels, n_cores)
%
% Modded by MDT -- temporary. Need to fix all of the parallel processing!
%
% IN PROGRESS!
%
% NOT FOR USE!


  %win_sz    = 10*dt; %.2        % seconds.
  n_win_sz  = 1; %floor( win_sz / dt );
  n         = size( data_in, 2 );
  n_win     = floor( n / n_win_sz );
  if( n_win == 0 )
    fprintf( '\n\nShould not happen.\n\n' );
    keyboard
    [ data_out, ref_est, nn_ref_est ] = robustCAR_helper( data_in, dt, n_min_contributing_channels, nm, taxis );
  else
    data_out    = zeros( size( data_in, 1 ), n );
    ref_est     = zeros( 1, n );
    nn_ref_est  = zeros( 1, n );
    fprintf( '\n\nrobustCAR4: Entering for jj = 1 : n ...\n\n' )
    tic
    if matlabpool( 'size' ) == 0
      matlabpool( n_cores );
    end
    parfor jj = 1 : n
      if mod( jj, 100 ) == 0 
        fprintf( 'Sample No: %d\n', jj )
      end
      i_start = jj; i_stop  = i_start + n_win_sz - 1;
      i_stop  = min( [ n i_stop ] );
      inds    = i_start:i_stop;
      if std( data_in(:,inds)) > 1e-20
        [ tmp_out, ref_out, nn_out ]  = robustCAR_helper( data_in(:,inds), dt, n_min_contributing_channels );

try
        data_out(:,jj)              = tmp_out(1,:)';   
        ref_est(jj)                 = ref_out(1); 
        nn_ref_est(jj)              = nn_out(1);
catch me
        me.message
        dbstack
        keyboard
end

      else
         data_out(:,jj)              = 0; ref_est(jj) = mean(data_in(:,jj)); nn_ref_est(jj) = size( data_in, 1 );
      end
    end % for jj 
    matlabpool close
    fprintf( '\n\nrobustCAR4_winsz_1b: Exiting for jj = 1 : n ...\n\n' )
    toc
    data_out = data_out';
  end
end

