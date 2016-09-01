function plot_ica_analysis(tsStack, conf, i_range, j_range)
%
% Here only tsStack is mandatory.

% SETTINGS
if (nargin < 2), conf.present = 1; end

% ensure standard settings are present at least
if (~isfield(conf, 'do_subtract_ref')),     conf.do_subtract_ref = 1; end
if (~isfield(conf, 'do_subtract_global')),  conf.do_subtract_global  = 1; end
if (~isfield(conf, 'do_color_one_only')),   conf.do_color_one_only   = 0; end
if (~isfield(conf, 'do_emphasize_col')),    conf.do_emphasize_col    = 0; end
if (~isfield(conf, 'ind_ref')),             conf.ind_ref             = 1:11; end
if (~isfield(conf, 'n_ica_comp')),          conf.n_ica_comp          = 3; end
if (~isfield(conf, 'w_bg')),                conf.w_bg                = 0.1; end
if (~isfield(conf, 'my_fontsize')),         conf.my_fontsize         = 14; end
if (~isfield(conf, 'lowres')),              conf.lowres              = 150; end
if (~isfield(conf, 'col')),                 conf.col                 = @hsv; end
%if (~isfield(conf, 'ind_skip')),            conf.ind_skip           = []; end
if (~isfield(conf, 'ind_skip_start')),      conf.ind_skip_start      = 1; end
if (~isfield(conf, 'ind_skip_end')),        conf.ind_skip_end        = 1; end

if (nargin < 3), i_range = 1:size(tsStack,1); end % If no argument, use the existing height dimensions
if (nargin < 4), j_range = 1:size(tsStack,2); end % If no argument, use the existing width dimensions

do_subtract_ref 	= conf.do_subtract_ref;
do_subtract_global  = conf.do_subtract_global;
do_color_one_only   = conf.do_color_one_only;
do_emphasize_col    = conf.do_emphasize_col;
ind_ref             = conf.ind_ref;
n_ica_comp          = conf.n_ica_comp;
w_bg                = conf.w_bg;
my_fontsize         = conf.my_fontsize;
lowres              = conf.lowres;
col                 = conf.f_col(n_ica_comp);
skipping_start      = conf.ind_skip_start;
skipping_end        = conf.ind_skip_end;

col = col / max(sum(col,1));
scale  = min(lowres / numel(i_range), lowres / numel(j_range));

% Convert from tiff to cropped and downsampled matrices
for c = 1:size(tsStack, 3)
    tmp = double(tsStack(i_range,j_range,c));
    tmp = imresize(tmp, scale, 'bilinear');
    tmp = tmp(2:end-1,2:end-1);
    
    if (c == 1)
        I = zeros(size(tmp,1), size(tmp,2), size(tsStack, 3));
    end
    
    I(:,:,c) = tmp;
end

% Pull out reference and mean
I_in   = I;
I_ref  = mean(I(:,:,ind_ref), 3);
I_mean = mean(I,3);
I_ref_max = max(I_ref(:));

S_noise = mean(mean(std(I(:,:,ind_ref),[],3)));

% Skip if requested
ind = ones(1,size(I,3));
ind(skipping_start) = 0;
ind(skipping_end) = 0;
I = I(:,:,ind == 1);
I_in = I_in(:,:,ind == 1);

if (do_subtract_ref)
    for c = 1:size(I,3)
        I(:,:,c) = I(:,:,c) - I_ref;
    end
end

if (do_subtract_global) % remove global signal mean
    for c = 1:size(I,3)
        I(:,:,c) = I(:,:,c) - mean(mean(I(:,:,c)));
    end
end

% Perform ICA analysis
data_in = reshape(I_in, size(I,1) * size(I,2), size(I, 3));
data    = reshape(I, size(I,1) * size(I,2), size(I, 3));

% center data to make projections useful
if (~do_subtract_ref)
    data = data - repmat(mean(data,2), 1, size(data,2));
end

[icasig, A, W] = fastica (data', ...
    'approach', 'symm', ...
    'g', 'pow3', ...
    'lastEig', n_ica_comp);

% Turn ICA components into normalized bases
Ac = A ./ repmat(sqrt(diag(A' * A))', size(A,1), 1);

% Flip sign of icasig to match data
for c = 1:size(A,2)
    tmp = data * Ac(:,c);
    tmp(abs(tmp) < 10 * S_noise) = 0;
    if (mean(tmp) < 0)
        A(:,c) = -A(:,c);
        Ac(:,c) = -Ac(:,c);
    end
end


% Order by elicited strength
mean_score = zeros(1,size(A,2));
for c = 1:size(A,2)
    tmp = data * Ac(:,c);
    mean_score(c) = mean(tmp);%(abs(tmp) > 10 * S_noise));
end

[~,ind] = sort(mean_score, 'descend');
A  = A(:,ind);
Ac = Ac(:,ind);
mean_score(ind)

% XXX: Cleaning here?
if (1)
    Ac(:, isnan(mean_score)) = 0;
    A(:, isnan(mean_score)) = 0;
end

% Calculate a score based on projection onto the normalized bases
qq = reshape(data * Ac, size(I,1), size(I, 2), size(A,2));
qq = abs(qq); % treat positive and negative projs equally
qq = qq - quantile(qq(:), 0.05);
qq = qq / quantile(qq(:), 0.9999);

if (do_emphasize_col) % suppress low values
    qq = (qq * 1.4).^2;
end

% Preparate a background for the display
I_bg = I_ref;
I_bg = I_bg - quantile(I_bg(:), 0.01);
I_bg = I_bg / quantile(I_bg(:), 0.99);

% Prepare a colorized images
p = zeros(size(I,1), size(I,2), 3);
for i = 1:1:size(I,1)
    for j = 1:1:size(I,2)
        w = squeeze(qq(i,j,:))';
        w(w < 0) = 0;
        
        % emphasize strongest component
        if (do_emphasize_col)
            w_max = max(w);
            w = w / w_max;
            w = w.^2;
            w = w * w_max;
        end
        
        if (do_color_one_only)
            w(w ~= max(w(:))) = 0;
            w = w.^(3/2);
        end
        
        p(i,j,:) = ...
            (0+w_bg) * I_bg(i,j) + ...
            (1-w_bg) * w * col;
    end
end

p(p > 1) = 1;
p(p < 0) = 0;

% Prepare time series plot
tmp = []; tmp_in = [];
n_plot = 8;
for c = 1:n_plot
    u = (c-1)/(n_plot-1) * size(data,2);
    x = 1:size(data,2);
    s = size(data,2) / n_plot * 1.5;
    w = exp(- (x - u).^2 / s^2);
    w = w / sum(w);
    tmp_in = [tmp_in reshape(data_in * w', size(I,1), size(I,2))];
    tmp = [tmp reshape(data * w', size(I,1), size(I,2))];
end

% scale for viewing
tmp = tmp / quantile(abs(tmp(:)), 0.99);
tmp = (tmp + 1) / 2;
tmp_in = tmp_in - quantile(abs(tmp_in(:)),0.01);
tmp_in = tmp_in / quantile(abs(tmp_in(:)), 0.99);

% Create a new window for the ICA
newICAfig=figure('Name', 'Independent component analysis','NumberTitle','off');

% replace print button with print preview instead, allowing for postscript export
hToolbar = findall(gcf,'tag','FigureToolBar');
hPrintButton = findall(hToolbar,'tag','Standard.PrintFigure');
set(hPrintButton, 'ClickedCallback','printpreview(gcbf)', 'TooltipString','Print Preview or Export as Postscript');

% plot time series
axes('position', [0.05 0.05 0.9 0.35]);
imagesc([tmp_in; tmp]);
axis image;
caxis([0 1]);
set(gca,'ytick', []);
set(gca,'xtick', [1:size(I,2):size(tmp,2) size(tmp,2)]);
set(gca,'xticklabel', num2str( round( (0:n_plot) / (n_plot ) * size(data,2))' ));

% plot colorized image
axes('position', [0.05 0.45 0.48 0.5]);
imagesc(p);
axis image off;

% plot curves
axes('position', [0.63 0.5 0.3 0.4]);
for c = 1:size(A,2)
    % If the Curve Fitting Toolbox is present it is possible to substitute the
    % fastsmooth function with smooth.
    %plot(smooth(A(:,c),3) / max(data(:)) * 100, 'k-', 'color', col(c,:) / max(col(:)) * 0.8, 'linewidth', 2); hold on;
    plot(fastsmooth(A(:,c),3) / max(data(:)) * 100, 'k-', 'color', col(c,:) / max(col(:)) * 0.8, 'linewidth', 2); hold on;
end
set(gca,'ytick',round([-1:0.1:1] * 100));
set(gca,'xtick',0:5:size(A,1));
box off;
xlabel('Frame', 'fontsize',  my_fontsize);
ylabel('Response in percent of max', 'fontsize',  my_fontsize);
legend(num2str((1:size(A,2))'));
perc_max = (floor(max(abs(A(:))) / max(data(:)) * 10) + 1) * 10;
axis([0 size(A,1)+1 [-1 1] * perc_max]);
set(gca,'fontsize', my_fontsize);