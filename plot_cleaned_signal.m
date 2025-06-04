function [clean_signal, clean_time] = plot_cleaned_signal(time, signal, sig_clean, bad_mask, chan_labels, asr_mode, vis_asr, max_fraction_bad)

% PLOT_CLEANED_SIGNAL: Visualize and optionally remove/reconstruct signal segments
%
% Inputs:
%   time        - Time vector [1 x N]
%   signal      - Signal matrix [channels x N]
%   sig_clean   - Reconstructed (cleaned) signal [channels x N], used only in 'reconstruct' mode
%   bad_mask    - Logical matrix [samples x channels]
%   chan_labels - Cell array of channel labels
%   asr_mode    - 'reconstruct' or 'remove'
%   vis_asr     - true/false (whether to plot)
%   max_fraction_bad - threshold for discarding sample if too many channels are bad
%
% Outputs:
%   clean_signal - Cleaned signal [channels x kept_samples] (only differs from input if remove mode)
%   clean_time   - Corresponding time vector (only if remove)

if nargin < 8
    max_fraction_bad = 0.5;
end

if size(bad_mask, 1) ~= length(time)
    bad_mask = bad_mask';  % ensure [samples x channels]
end

signal = double(signal);  % enforce numeric
[num_chans, num_samples] = size(signal);
clean_signal = signal;
clean_time = time;

if vis_asr
    % Normalize for consistent vertical range
    signal = (signal - mean(signal, 2, 'omitnan')) ./ std(signal, 0, 2, 'omitnan');
    if strcmpi(asr_mode, 'reconstruct') && ~isempty(sig_clean)
        sig_clean = (sig_clean - mean(sig_clean, 2, 'omitnan')) ./ std(sig_clean, 0, 2, 'omitnan');
    end
    offset_step = 3;
end

if strcmpi(asr_mode, 'reconstruct')
    clean_signal = sig_clean;
    clean_time = time;

    if vis_asr
        figure('color','w'); hold on

        for iChan = 1:num_chans
            offset = (num_chans - iChan) * offset_step;
            t = time(:)';

            % Plot reconstructed signal (blue, always visible)
            scrollplot({t, sig_clean(iChan,:) + offset, 'color', '#0072BD'}, {'X'}, 15);

            % Plot original signal only at bad samples (red)
            bad_idx = find(bad_mask(:, iChan));
            if ~isempty(bad_idx)
                diffs = diff(bad_idx) > 1;
                breaks = [1, find(diffs)+1, numel(bad_idx)+1];
                for s = 1:(length(breaks)-1)
                    seg = bad_idx(breaks(s):breaks(s+1)-1);
                    seg_x = [t(seg), NaN];
                    seg_y = [signal(iChan, seg) + offset, NaN];
                    scrollplot({seg_x, seg_y, 'color', 'r'}, {'X'}, 15);
                end
            end
        end

        yticks((0:num_chans-1) * offset_step);
        yticklabels(flip(chan_labels(:)));
        ylim([-1 offset_step * num_chans]);
        title('Signal Cleaning Overview');
    end


elseif strcmpi(asr_mode, 'remove')
    keep_idx = sum(bad_mask, 2) <= max_fraction_bad * size(bad_mask, 2);

    if vis_asr
        figure('color','w'); hold on
        bad_all_x = {}; bad_all_y = {};
        clean_all_x = {}; clean_all_y = {};

        for iChan = 1:num_chans
            signal_i = signal(iChan, :);
            time_i = time(:)';
            offset = (num_chans - iChan) * offset_step;

            bad_idx = find(bad_mask(:, iChan));
            good_idx = ~bad_mask(:, iChan);

            trace = nan(1, length(signal_i));
            trace(good_idx) = signal_i(good_idx) + offset;
            clean_all_x{end+1} = time_i;
            clean_all_y{end+1} = trace;

            if ~isempty(bad_idx)
                diffs = diff(bad_idx) > 1;
                breaks = [true, diffs(:).', true];
                edges = find(breaks);
                for s = 1:(length(edges) - 1)
                    segment = bad_idx(edges(s):(edges(s+1)-1));
                    bad_all_x{end+1} = [time_i(segment), NaN];
                    bad_all_y{end+1} = [signal_i(segment) + offset, NaN];
                end
            end
        end
        
        scrollplot({horzcat(clean_all_x{:})', horzcat(clean_all_y{:})', 'color', '#0072BD'}, {'X'}, 15);
        if ~isempty(bad_all_x)
            scrollplot({horzcat(bad_all_x{:})', horzcat(bad_all_y{:})', 'color', 'r'}, {'X'}, 15);
        end
        
        yticks((0:num_chans-1) * offset_step);
        yticklabels(flip(chan_labels(:)));
        ylim([-1 offset_step * num_chans]);
        title('Signal Cleaning Overview');
    end

    clean_signal = signal(:, keep_idx');
    clean_time = time(keep_idx');
else
    error('Unsupported asr_mode. Use ''reconstruct'' or ''remove''.');
end
