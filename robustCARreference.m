function [data_out, ref_est, nn_ref_est] = robustCARreference(data_in, dt, n_min_contributing_channels, n_cores)

% [data_out ref_est nn_ref_est] = robustCARreference(data_in, dt, n_min_contributing_channels, n_cores)
%
% This function is a modified and fixed version of the functions developed
% and released onto the web by Kyle Q. Lepage, Boston University. Obtained
% from the link: http://math.bu.edu/people/lepage/code.html during 9/2016.
% The downloaded code did not run (functions were not named the same as
% their filenames and there were other issues. Hopefully this version fixes
% these minor bugs and only introduces some code and efficiency
% optimizations. Hopefully. Neither the current editor (Matthew Turner) nor
% the original author accept any responsibility for anything having to do
% with this code or its use for *anything*. Please see original author's
% remarks at the bottom of this file.
%
% This code is unlicensed as the original author did not assign one. All
% original code in this repository is MIT licensed. Please see the license
% file for this project for details.
%
% Mandatory parameters:
%
% data_in   - A channels X time (samples) matrix
% dt        - The time length of a sample (1/Fs or sr)
%
% Optional parameters:
%
% n_min_contributing_channels - Minimum channels for calculation (8)
% n_cores                     - Number of cores on processor to use (4)
%
% Returns the rCAR re-referenced signals.

    if nargin < 3
        n_min_contributing_channels = 8;    % For us, 8 out of 14
        n_cores                     = 4;    % For testing, should be 7
    end
    
    % Initialization

    [n_channels, n_times] = size(data_in);
    z                     = 2^ceil(log2(n_times));
    half_z                = z/2;
    t                     = [0:n_times - 1]' * dt;
    
    % Remove the sample average.
    
    mean_data_in = mean(data_in, 2);
    data_in      = data_in - mean_data_in * ones(1,n_times);
    
    % ==============================================================
    % Perform in frequency domain as opposed to time domain.
    % ==============================================================
    data                      = data_in';
    fdata                     = fft( data, z );
    re_fdata                  = real( fdata( 1 : half_z+1, : ))';
    im_fdata                  = imag( fdata( 1 : half_z+1, : ))';
    
    [~, re_ref_est, nn_re_ref_est]     = robustCAR4_winsz( re_fdata, dt, n_min_contributing_channels, n_cores );
    [~, im_ref_est, nn_im_ref_est]     = robustCAR4_winsz( im_fdata, dt, n_min_contributing_channels, n_cores );
    
    re_ref_est_full           = [ re_ref_est fliplr( re_ref_est( 2 : end - 1 )) ];
    im_ref_est_full           = [ im_ref_est -fliplr( im_ref_est( 2 : end - 1 )) ];
    ref_est_full              = re_ref_est_full + 1i * im_ref_est_full;
    ref_est                   = ifft( ref_est_full );
    ref_est                   = ref_est( 1 : n_times );
    data_out                  = data_in - ones(n_channels,1) * ref_est;
    data_out                  = data_out';
    nn_ref_est                = ( nn_re_ref_est + nn_im_ref_est ) / 2;
    
end

% Original authors comments:
%
% Author:     Kyle Q. Lepage
%             The author takes no responsibility for anything that results
%             from the use of this code for anything.
%
% data_in -  n_channels  x n_times
% dt      - sample period in seconds.
%
% n_min_contributing_channels - minimum # of channels to contribute to
%                               reference estimate. 
% n_cores                     - number of parallel processes to run.
%
% Side-effects:
%             - removes the sample average across time.