% data_in     -  n_channels  x n_times
%
% 
% 
% Side Effect:    Removes the mean of the data.
%
function [ data_out ref_est nn_ref_est ] = ...
            robustCAR3c( data_in, dt, n_min_contributing_channels, nm, taxis )

  % ===================================================
  % Set the bisquare function cutoff, c
  % ===================================================
  c         = 2.0;

  if nargin < 3 
    n_min_contributing_channels = 8;
  end

  b_mad_scale_est         = true;
if exist( 'nm' )                                            
  fprintf( '\n\n============================================================\n' );
  if( b_mad_scale_est )
    fprintf( 'b_mad_scale_est set to true.\n' ); 
  else
    fprintf( 'b_mad_scale_est set to false.\n' ); 
  end
  fprintf( '\n============================================================\n' );
end

  [ n_channels n_times ]    = size( data_in );
  t                         = [ 0 : n_times - 1 ]' * dt;

if( ~b_mad_scale_est )
  % ==============================================
  % Robust MAR.
  % ==============================================
  data_in   = data_in - ( mean( data_in, 2 ) * ones(1,n_times) );
  %figure(1);clf;
  %  imagesc(data_in), colorbar,title('centered data')
  %  print -depsc2 ../out/d.eps, close(1);
  ts                      = [];
  ts.bDebug               = 0;
  ts.MW.order             = 1;
  ts.MAR.order            = 2 * ts.MW.order + 1;
  ts.dTs                  = dt;
  ts.MAR.bUseBothSides    = 1;
  ts.MAR.bUseCentreVals   = 0;
  for j = 1 : n_channels
    ts.data{j}  = data_in(j,:)';
    %figure(1);clf;plot( ts.data{j}, '-k' ), print( 'out/ts_d.ps',  '-dpsc2',  '-append' ); close(1);
  end
  ts = tsEstRobustMAR( ts );

  ts.MAR.order  = ( ts.MAR.order - 1 ) / 2;
  try 
    ts            = tsEstMultivarCovPastFuture( ts );
  catch me
    me.message
    dbstack
    keyboard
  end
  ts.MAR.scale  = median( abs( ts.MAR.d - ts.MAR.bothsides_est ));

  %nm = 'out/mar.eps';
  if exist( 'nm' )
  figure(1);clf;
    subplot(2,2,1),plot( ts.MAR.scale, '+k' ); title( 'scale' )
    subplot(2,2,3),imagesc( data_in ); colorbar, title( 'data' )
    subplot(2,2,2),imagesc( ts.MAR.bothsides_est' ); colorbar, title( 'bothsides MAR estimate' )
    subplot(2,2,4),imagesc( data_in - ts.MAR.bothsides_est' ); colorbar, title( 'bothside fit: residuals' )
    print( nm,  '-dpsc2', '-append' ), close(1);
  end


  % Clean the data with the MAR.
  ts            = tsCleanDataWithRobustMAR( ts );

  %clean_reshape_day( r, : ) = rccTS.MAR.cleaned_d( :, 1 )';
  ts.MAR.order  = 2 * ts.MAR.order + 1;
else
  if n_times == 1 
    abs_dev               = abs( data_in - median( data_in ));
    abs_dev               = abs_dev( 1:end-1 );
    ts.MAR.scale          = mean( abs_dev ) * ones(1,n_channels);
    if any( ~isfinite( ts.MAR.scale ))
      fprintf( '\n\nrobustCAR3b(): NaN detected.  Check scale?\n\n' );
      dbstack
      keyboard
    end
  elseif n_times < 15
    abs_dev               = abs( data_in - ( median( data_in, 2 )*ones(1,n_times) ) );

    % Trim the high value.
    abs_dev               = sort( abs_dev, 2 );
    abs_dev               = abs_dev( :, 1:end-1 );
    ts.MAR.scale          = mean( abs_dev, 2 )';
    i_bad                 = find( ts.MAR.scale < 1e-12 );
    ts.MAR.scale( i_bad ) = 1e-12;
  else
    ts.MAR.scale  = median( abs( data_in - ( median( data_in, 2 )*ones(1,n_times) ) ), 2 )';   % These are often much larger than the 
  end
end


  % ==================================================================================
  % Robust M-estimate of the mean using a scale =  ts.MAR.scale  
  % (estimated innovations scale)
  % ==================================================================================
  %ts.MAR.scale = ts.MAR.scale;
  n_max_iters   = 10;
  sse_threshold = 1e-30;
  n_starts      = 1; %10;
  mu            = zeros( n_times, n_max_iters, n_starts );
  nnz           = zeros( n_times, n_max_iters, n_starts );
  sse           = 1e10*ones( n_times, n_max_iters, n_starts );
  ref_est       = zeros( n_times, 1 );
  for ii = 1 : n_times
    d_this_time       = data_in( :, ii );
    median_this_time  = median( d_this_time );
    for i_start = 1 : n_starts 
      %mu( ii, 1, i_start ) = .1 * ts.MAR.scale(1) * randn(1) + median_this_time;
      mu( ii, 1, i_start )    = median_this_time;
      x                       = ( d_this_time - mu( ii, 1, i_start ) ) ./ ts.MAR.scale';
      bx                      = bisquarePsiFunction( x, c );
      n_non_zero              = length( bx( abs( bx ) > 1e-14 ));
      if( n_non_zero > 0 )
        e( ii, 1, i_start )     = sum( bisquarePsiFunction( x, c ));
        sse( ii, 1, i_start )   = e( ii, 1, i_start )^2 / n_non_zero; % ~ chisquared 1 under null
      else
        e( ii, 1, i_start )     = 1e10; %sum( bisquarePsiFunction( x ));
        sse( ii, 1, i_start )   = 1e10; %e( ii, 1, i_start )^2 / n_non_zero; % ~ chisquared 1 under null
      end

      kk = 2;
      while sse( ii, kk-1, i_start ) > sse_threshold || ~isfinite( sse( ii, kk-1, i_start )) 
        if kk > n_max_iters
          fprintf( '\n\nMax iters reached. Breaking.  Min( SSE ) = %.3e\n\n', min( sse(:) ) );
%          keyboard
          break;
        end
        base_delta              = get_delta( d_this_time, mu( ii, kk-1, i_start ), ts.MAR.scale', c );
        try_this_mu             = mu( ii, kk-1, i_start ) + base_delta;
        x                       = ( d_this_time - try_this_mu ) ./ ts.MAR.scale';
        n_non_zero              = sum( abs( bisquarePsiFunction( x, c )) > 0 );
        mu( ii, kk, i_start )   = try_this_mu;
        e( ii, kk, i_start )    = sum( bisquarePsiFunction( x, c ));
        nnz( ii, kk, i_start )  = n_non_zero;
        sse( ii, kk, i_start )  = e( ii, kk, i_start )^2 / n_non_zero; % ~sum of squared standard normal.
        kk                      = kk + 1;
      end  % while
    end % istart

    % Assumed: only 1 start.
    i_use                     = find( isfinite( sse( ii, :, 1 )));
    [ min_sse  i_min ]        = min( sse( ii, i_use, 1 ) );
    i_iter                    = i_use( i_min );

    if( min_sse > 1e-4 ) %|| isnan( min_end_time_abs_e_over_starts ))
      fprintf( '\n\nnot converged\n\n' );
      dbstack
      keyboard
    else
      nn_ref_est( ii )                                  = nnz( ii, i_iter, i_start );
      ref_est( ii )                                     = mu( ii, i_iter, i_start );
      err( ii )                                         = sse( ii, i_iter, i_start );
      used_iter( ii )                                   = i_iter;
      if isnan( err(ii) )
        keyboard
        ref_est( ii ) = median( d_this_time );
        x             = ( d_this_time - ref_est(ii) ) ./ ts.MAR.scale';
        err( ii )     = sum( bisquarePsiFunction( x, c ));
      end
      %fp = fopen( '../out/debug.txt', 'a+' );
      %fprintf( fp, 'Used iter = %d\tsse = %.2e\tn non zero = %d\n', i_iter, err(ii), nn_ref_est(ii) );
      %fprintf( 'Used iter = %d\tsse = %.2e\tn non zero = %d\n', i_iter, err(ii), nn_ref_est(ii) );
      %fclose( fp );
      if err( ii ) > 1e-4
        keyboard
      end
    end
  end % time index


  %x = -2 : .01 : 2;
  %figure(1);clf;
  %  plot( x, bisquarePsiFunction( x ), '-b' )  ; title( 'bisquare psi' ), hold on
  %  plot( x, deriv_bisquarePsiFunction( x ), '-g'   ); title( 'deriv. bisquare psi' ); hold on
  %  print -depsc2 ../out/chk.eps, close(1);

  % ==============================================
  % The reference estimate can have outliers.
  % Give it a cleaning.
  % ==============================================
  if( 0 )
    ts.data{1}      = ref_est;
    ts.AR.order     = 2;
    ts.AR.nIters    = 1;
    %[ ar2Est residuals ARCoeffs ] = PreWhitenData( ref_est, ts.AR.order, false )
    how_robust  = 'a little';
    ts          = tsEstRobustAR( ts, how_robust );
    ts          = tsCleanDataWithRobustAR( ts );
    ref_est     = ts.clean;
  end


  % ==============================================
  % Referenced data.
  % ==============================================
  data_out    = ( data_in - ones( n_channels, 1 ) * ref_est' ) ;

  if exist( 'nm' )
  figure(1);clf;
    subplot(3,2,1),plot( taxis, ref_est, '-k' ); title( 'reference signal estimate' );
    subplot(3,2,2),plot( taxis, err, '-k+' ); , title( [{'robust estimation error (squared & normalized)'},{''}] );
      xlabel( 'time' ), ylabel( 'error' )
    subplot(3,2,4),imagesc( data_out ); c =  caxis;
    colorbar, title( 'rCARIII data' )
    subplot(3,2,3),imagesc( taxis, 1:n_channels, data_in ); colorbar, title( 'data' ); caxis( c );
    subplot(3,2,5), plot( taxis, nn_ref_est, 'k+' ), colorbar, title( 'Number of Used Channels' ),  set(gca,'ylim',[-1 20] )
    subplot(3,2,6), plot( taxis, used_iter, 'k+' ), colorbar, title( 'Used Iteration' ), %set(gca,'xlim',[1 120] )
    %print -depsc2 ../out/rCARII.eps, close(1);
    print( nm,  '-dpsc2', '-append' ), close(1);
  end

  if( 0 )
  keyboard
  figure(2);clf;
  subplot(1,2,1), plot( taxis, median( data_in ), '-k' );title('median datain')
  subplot(1,2,2), imagesc( data_in - ( ones( 19, 1 ) * median( data_in ) ) ), title( 'median re-ref data' )
    print( nm,  '-dpsc2', '-append' ), close(2);
    end % if( 0 )

  data_out    = data_out';
end

function rvals = bisquarePsiFunction( arg, c )
     zinds           = find( abs( arg ) > c );
     rvals           = arg .* ( 1.0 - ( arg / c ).^2 ).^2;
     rvals( zinds )  = 0.0;
end

function rvals = deriv_bisquarePsiFunction( x, c )
   zinds           = find( abs( x ) > c );
   rvals           = 1 - 6 * (x/c).^2 + 5 * ( x/c ).^4;
   rvals( zinds )  = 0.0;

end

function r = dBPF_dmu( x, c )
  r           = deriv_bisquarePsiFunction( x,c );
end

function r = d2BPF_dmu2( x, c )
  r           = -12/c^2 * x + 20/c^4 * x.^3;
  zinds       = find( abs( x ) > c );
  r( zinds )  = 0;
end

function r = d3BPF_dmu3( x, c )
  r           = 60 * x.^2 / c^4 - 12 / c^2;
  zinds       = find( abs( x ) > c );
  r( zinds )  = 0;
end

function r = d4BPF_dmu4( x, c )
  r           = 120 * x / c^4;
  zinds       = find( abs( x ) > c );
  r( zinds )  = 0;
end

function r = d5BPF_dmu5( x, c )
  r           = 120 * ones(size(x)) / c^4;
  zinds       = find( abs( x ) > c );
  r( zinds )  = 0;
end


function   delta  = get_delta( d, mu, scales, c )

  if( 0 )
     d         = median( scales(:,1) ) * randn( 20, 1 );
    dmu       = .01 * median( scales(:));
    mu        = ones( size(d,1), 1 ) * [ -3*median(scales) : dmu : 3*median(scales) ];
    scales    = scales * ones(1, size(mu, 2));
    d         = d * ones( 1, size( mu, 2 ));
  end

  x           = ( d - mu ) ./ scales;
  n_non_z     = sum( bisquarePsiFunction(x,c) ~= 0 );
  e           = sum( bisquarePsiFunction( x,c ));
  de_dmu      = sum( dBPF_dmu( x,c ) ./ -scales );
  d2e_dmu2    = sum( d2BPF_dmu2( x,c ) ./ ( -scales ).^2 );
  d3e_dmu3    = sum( d3BPF_dmu3( x,c ) ./ ( -scales ).^3 );
  d4e_dmu4    = sum( d4BPF_dmu4( x,c ) ./ ( -scales ).^4 );
  d5e_dmu5    = sum( d5BPF_dmu5( x,c ) ./ ( -scales ).^5 );

  if( 0 )
    nvals     = size( mu, 2 )
    mux       = linspace( -3, 3, nvals );
    dx        = mux(2)-mux(1);


    figure(1);clf;
    plot(  mux, e, '-k+' );
    print -depsc2 out/chk.eps, close(1);

    figure(1);clf;
    plot(  mux(1:end-1), diff(e)/dmu , '-k+', mux, de_dmu, '-r' );
    %set(gca,'xlim', [-.01 .01] );
    print -depsc2 out/chk.eps, close(1);


    figure(1);clf;
    plot(  mux(1:end-1), diff(de_dmu)/dmu , '-k+', mux, d2e_dmu2, '-r' );
    print -depsc2 out/chk.eps, close(1);
    keyboard
    figure(1);clf;
    plot(  mux(1:end-1), diff( d2e_dmu2)/dmu , '-k+', mux, d3e_dmu3, '-r' );
    print -depsc2 out/chk.eps, close(1);
    figure(1);clf;
    plot(  mux(1:end-1), diff( d3e_dmu3)/dmu , '-k+', mux, d4e_dmu4, '-r' );
    print -depsc2 out/chk.eps, close(1);
    figure(1);clf;
    plot(  mux(1:end-1), diff( d4e_dmu4 )/dmu , '-k+', mux, d5e_dmu5, '-r' );
    print -depsc2 out/chk.eps, close(1);
  end

  % Prepare to call roots.
  c0          = 1/factorial(5) * d5e_dmu5;          c3        = 1/factorial(2) * d2e_dmu2;
  c1          = 1/factorial(4) * d4e_dmu4;          c4        = 1/factorial(1) * de_dmu;
  c2          = 1/factorial(3) * d3e_dmu3;          c5        = e;
  c_vec       = [ c0 c1 c2 c3 c4 c5 ];

  % Second derivative is discontinuous at 1.
  %c_vec       = [ c2 c3 c4 c5 ];
  if any( ~isfinite( c_vec ))
    fprintf( '\n\nrobustCAR3b(): NaN detected.  Check scale?\n\n' );
    dbstack
      keyboard
  end

  candidate_deltas                                  = roots( c_vec );
  candidate_deltas( imag( candidate_deltas ) ~= 0 ) = [];
  if ~isempty( candidate_deltas )

  try
    n_candidates      = length( candidate_deltas );
    new_mus           = mu + make_col( candidate_deltas )';
    x                 = ( d*ones(1,n_candidates) - ones(size(d,1),1)*new_mus ) ./ ( scales * ones(1,n_candidates));
    newe_s            = sum( bisquarePsiFunction( x,c ));
    [ tmp i_mu ]      = min( abs( newe_s ));
    newe              = newe_s( i_mu );
    mu                = new_mus( i_mu );
    delta_e           = abs( newe ) -abs( e);
    delta             = candidate_deltas( i_mu );
  catch me
    me.message
    keyboard
  end
    %fprintf( '\nDelta Abs err = %.3e   New abs error = %.3e\n\n', delta_e, abs( newe  ));
    n_non_z2          = sum( bisquarePsiFunction(x,c) ~= 0 );
    delta_n_non_z     = n_non_z - n_non_z2;
  else
    candidate_deltas  = 0;
    delta             = 0;
    new_mu            = median(d) + .2 * scales(1) * randn(1);
    delta             = new_mu - mu;
    fprintf( '\nAll candidates complex valued: mu-median%.2e\n', mu - median(d));
  end
end                                                                    


