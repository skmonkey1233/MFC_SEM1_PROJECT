%% ========================================================================
%  MINI PROJECT: Image Compression Using Singular Value Decomposition (SVD)
%  Course  : Mathematics for Computing (MRIIRS)
%  Units   : 1 (Matrices, rank), 2 (Eigenvalues, SVD), 5 (Optimization)
%  Software: MATLAB R2022a or later + Image Processing Toolbox
%
%  PARTS
%    A. Core project (grayscale image)         -> Sections 1-8
%    B. Extra 1: Colour (RGB) compression      -> Section 9
%    C. Extra 2: Several test images compared  -> Section 10
%    D. Extra 3: Interactive slider demo       -> Section 11
%  Run the whole file (F5). Close the slider window to finish.
%  ========================================================================
clc; clear; close all;

%% 1. LOAD IMAGE AND CONVERT TO GRAYSCALE ---------------------------------
imgFile = 'cameraman.tif';          % Image Processing Toolbox sample image
I = readImageOrSelect(imgFile, 'Select the image to compress');
if size(I,3) == 3
    I = rgb2gray(I);
end

%% 2. REPRESENT THE IMAGE AS A MATRIX -------------------------------------
A = double(I);
[m, n] = size(A);
if min(m, n) < 2
    error('svd_image_compression:ImageTooSmall', ...
          'Select an image that is at least 2 pixels high and 2 pixels wide.');
end

%% 3. APPLY SVD:  A = U * S * V' ------------------------------------------
[U, S, V] = svd(A, 'econ');          % reduced SVD: less memory/work
sigma = diag(S);
p = numel(sigma);

% Unit 1: rank = number of non-zero singular values (same tolerance as
% MATLAB's rank(), but reuses the SVD instead of computing a second one).
r = sum(sigma > max(m,n) * eps(sigma(1)));
fprintf('Image size          : %d x %d\n', m, n);
fprintf('Rank of image matrix: %d\n', r);

% Unit 2: singular values are square roots of the eigenvalues of A'*A.
% This independent check is useful for the report, but a second matrix
% factorisation is unnecessarily slow for a large, user-selected photo.
if max(m,n) <= 1024
    ev = sort(eig(A'*A), 'descend');
    fprintf('Check: max |sigma - sqrt(eig(A''A))| = %.3e\n', ...
            max(abs(sigma - sqrt(max(ev(1:p),0)))));
else
    fprintf('Eigenvalue check skipped for image dimensions above 1024 pixels.\n');
end

%% 4. RECONSTRUCT USING k SINGULAR VALUES (rank-k approximation) ----------
%  A_k = sum_{i=1..k} sigma_i * u_i * v_i'
kList = [5 10 20 50 100];
kList = kList(kList <= p);

kAll = 1:p;

% The discarded singular-value energy gives the error exactly, so there is
% no need to reconstruct a full image once for every possible value of k.
% This avoids p large matrix multiplications and makes this section fast.
singularEnergy = sigma.^2;
totalEnergy    = sum(singularEnergy);
energy         = cumsum(singularEnergy) / totalEnergy;
tailEnergy     = flipud(cumsum(flipud(singularEnergy)));
discardedEnergy = [tailEnergy(2:end); 0];
relErr         = sqrt(discardedEnergy / totalEnergy).';
mseVal         = (discardedEnergy / (m*n)).';
psnrVal        = 10*log10(255^2 ./ max(mseVal, eps));

%% 5. STORAGE REQUIREMENTS AND COMPRESSION --------------------------------
%  Original : m*n numbers
%  Rank-k   : k*(m + n + 1) numbers
storageOrig = m*n;
storageK    = kAll * (m + n + 1);
CR          = storageOrig ./ storageK;
savingPct   = (1 - storageK/storageOrig) * 100;
kBreakEven  = floor(m*n / (m + n + 1));

fprintf('\nBreak-even rank (k <= %d gives real compression)\n\n', kBreakEven);
fprintf('%6s %12s %10s %10s %10s %10s\n', ...
        'k', 'Storage', 'CR', 'Saving%', 'RelErr', 'PSNR(dB)');
for k = kList
    fprintf('%6d %12d %10.2f %10.2f %10.4f %10.2f\n', ...
        k, storageK(k), CR(k), savingPct(k), relErr(k), psnrVal(k));
end

k99   = find(energy >= 0.99, 1);
k30dB = find(psnrVal >= 30, 1);
fprintf('\nk for 99%% energy      : %d (CR = %.2f)\n', k99, CR(k99));
fprintf('k for PSNR >= 30 dB   : %d (CR = %.2f)\n', k30dB, CR(k30dB));

%% 6. DISPLAY ORIGINAL AND COMPRESSED IMAGES ------------------------------
figure('Name','Original vs Compressed','Position',[50 50 1100 650]);
nPlots = numel(kList) + 1;
subplot(2, ceil(nPlots/2), 1);
imshow(uint8(A)); title(sprintf('Original (rank %d)', r));
for i = 1:numel(kList)
    k  = kList(i);
    Ak = reconstructRankK(U, sigma, V, k);
    subplot(2, ceil(nPlots/2), i+1);
    imshow(uint8(min(max(Ak,0),255)));
    title(sprintf('k = %d, CR = %.1f:1, PSNR = %.1f dB', ...
          k, CR(k), psnrVal(k)), 'FontSize', 9);
end

%% 7. GRAPHS: SINGULAR VALUES, ERROR, ENERGY, TRADE-OFF -------------------
figure('Name','Error Analysis','Position',[80 80 1100 700]);

subplot(2,2,1);
semilogy(kAll, sigma, 'b-', 'LineWidth', 1.5); grid on;
xlabel('Index i'); ylabel('\sigma_i (log scale)');
title('Singular value spectrum');

subplot(2,2,2);
plot(kAll, relErr, 'r-', 'LineWidth', 1.5); grid on;
xlabel('Number of singular values k');
ylabel('||A-A_k||_F / ||A||_F');
title('Reconstruction error vs k');

subplot(2,2,3);
plot(kAll, energy*100, 'g-', 'LineWidth', 1.5); hold on;
yline(99, 'k--', '99% energy'); grid on;
xlabel('Number of singular values k'); ylabel('Cumulative energy (%)');
title('Energy retained vs k');

subplot(2,2,4);
plot(savingPct, psnrVal, 'm-', 'LineWidth', 1.5); grid on;
xlabel('Storage saved (%)'); ylabel('PSNR (dB)');
title('Quality vs storage saving (trade-off)');
xlim([0 100]);

%% 8. OPTIMIZATION VIEW (Unit 5) ------------------------------------------
%  Eckart-Young theorem: A_k solves
%      minimize ||A - B||_F   subject to  rank(B) <= k.
%  Minimum error = sqrt(sum of discarded sigma_i^2).
k = min(50, p-1);
theory = sqrt(sum(sigma(k+1:end).^2));
Ak = reconstructRankK(U, sigma, V, k);
fprintf('\nEckart-Young check (k=%d): ||A-Ak||_F = %.4f, theory = %.4f\n', ...
        k, norm(A-Ak,'fro'), theory);

targetErr = 0.05;                   % allowed relative error (the constraint)
kOpt = find(relErr <= targetErr, 1);
fprintf('Smallest k with relative error <= %.0f%% : k = %d (CR = %.2f)\n', ...
        targetErr*100, kOpt, CR(kOpt));

%% 9. EXTRA 1: COLOUR (RGB) IMAGE COMPRESSION -----------------------------
%  Apply SVD separately to the R, G and B channel matrices.
C = double(readImageOrSelect('peppers.png', ...
           'Select an RGB image for the colour-compression demo'));
if size(C,3) == 1
    C = repmat(C, 1, 1, 3);           % keep the RGB demo valid for grayscale input
elseif size(C,3) > 3
    C = C(:,:,1:3);                    % discard an alpha channel, if present
end
[mc, nc, ~] = size(C);
kc = [5 20 50];
kc = kc(kc <= min(mc, nc));
[Uc, Sc, Vc] = deal(cell(1,3));
sc = cell(1,3);
for ch = 1:3
    [Uc{ch}, Sc{ch}, Vc{ch}] = svd(C(:,:,ch), 'econ');
    sc{ch} = diag(Sc{ch});
end
figure('Name','Colour compression','Position',[60 60 1100 330]);
subplot(1, numel(kc)+1, 1);
imshow(uint8(C)); title('Original RGB');
for i = 1:numel(kc)
    k  = kc(i);
    Ck = zeros(size(C));
    for ch = 1:3
        Ck(:,:,ch) = reconstructRankK(Uc{ch}, sc{ch}, Vc{ch}, k);
    end
    Ck   = min(max(Ck,0),255);
    mseC = mean((C(:) - Ck(:)).^2);
    crC  = (3*mc*nc) / (3*k*(mc + nc + 1));   % 3 channels cancel out
    subplot(1, numel(kc)+1, i+1);
    imshow(uint8(Ck));
    title(sprintf('k = %d, CR = %.1f:1, PSNR = %.1f dB', ...
          k, crC, 10*log10(255^2/max(mseC,eps))), 'FontSize', 9);
end

%% 10. EXTRA 2: COMPARE SEVERAL TEST IMAGES -------------------------------
%  Smooth images have fast-decaying singular values (compress well);
%  detailed/textured images have slow decay (compress poorly).
files = {'cameraman.tif', 'circuit.tif', 'text.png'};
kTest = 20;
figure('Name','Image comparison','Position',[70 70 1100 700]);
subplot(2,1,1); hold on; grid on;
fprintf('\n%-16s %8s %12s %12s %10s\n', ...
        'Image', 'Size', 'k(99% E)', 'RelErr', 'PSNR (dB)');
for f = 1:numel(files)
    T = readImageOrSelect(files{f}, ...
        sprintf('Select an image to use instead of %s', files{f}));
    if size(T,3) == 3, T = rgb2gray(T); end
    T = double(T);
    [mt, nt] = size(T);
    [Ut, St, Vt] = svd(T, 'econ');
    st = diag(St);
    semilogy(st/st(1), 'LineWidth', 1.5);       % normalised spectrum
    et   = cumsum(st.^2)/sum(st.^2);
    k99t = find(et >= 0.99, 1);
    kt   = min(kTest, numel(st));
    Tk   = reconstructRankK(Ut, st, Vt, kt);
    ee   = norm(T-Tk,'fro')/norm(T,'fro');
    pp   = 10*log10(255^2/mean((T(:)-Tk(:)).^2));
    fprintf('%-16s %4dx%-4d %12d %12.4f %10.2f\n', ...
            files{f}, mt, nt, k99t, ee, pp);
end
set(gca,'YScale','log');
xlabel('Index i'); ylabel('\sigma_i / \sigma_1 (log scale)');
title('Normalised singular value decay for different images');
legend(files, 'Interpreter','none');
subplot(2,1,2); axis off;
text(0.02, 0.5, ['Slower decay = more detail = needs larger k for the same ' ...
     'quality (lower compression).'], 'FontSize', 12);

%% 11. EXTRA 3: INTERACTIVE SLIDER ----------------------------------------
%  Move the slider to change k and watch image, error and storage update.
fig = figure('Name','Interactive SVD compression','Position',[100 100 800 620]);
ax  = axes('Parent', fig, 'Position', [0.1 0.2 0.8 0.72]);
imshow(uint8(A), 'Parent', ax);
txt = uicontrol('Style','text','Units','normalized', ...
      'Position',[0.1 0.10 0.8 0.05], 'FontSize', 11);
sld = uicontrol('Style','slider','Units','normalized', ...
      'Position',[0.1 0.04 0.8 0.04], 'Min',1, 'Max',p, 'Value',min(20,p), ...
      'SliderStep',[1/(p-1) 10/(p-1)]);
sld.Callback = @(src,~) updateView(src, ax, txt, U, sigma, V, A, ...
                                   relErr, psnrVal, CR);
updateView(sld, ax, txt, U, sigma, V, A, relErr, psnrVal, CR);

%% ---------------------------- LOCAL FUNCTIONS ---------------------------
function updateView(src, ax, txt, U, sigma, V, A, relErr, psnrVal, CR)
    k  = max(1, round(src.Value));
    Ak = reconstructRankK(U, sigma, V, k);
    imshow(uint8(min(max(Ak,0),255)), 'Parent', ax);
    txt.String = sprintf('k = %d | CR = %.2f:1 | Rel. error = %.4f | PSNR = %.1f dB', ...
                         k, CR(k), relErr(k), psnrVal(k));
end

function Ak = reconstructRankK(U, sigma, V, k)
% Form U(:,1:k)*diag(sigma(1:k))*V(:,1:k)' without constructing diag(...).
    Ak = (U(:,1:k) .* sigma(1:k).') * V(:,1:k)';
end

function I = readImageOrSelect(defaultFile, dialogTitle)
% Read an image and return it as uint8 (grayscale or RGB, values 0..255).
% Search order: 1) images/ folder next to this script,
%               2) MATLAB sample image on the path,
%               3) ask the user to select a local image.
    here = fileparts(mfilename('fullpath'));
    if isempty(here), here = pwd; end
    localCopy = fullfile(here, 'images', defaultFile);

    if isfile(localCopy)
        file = localCopy;
    elseif exist(defaultFile, 'file') == 2
        file = defaultFile;
    else
        file = selectImageFile(defaultFile, dialogTitle);
    end

    [I, map] = imread(file);
    if ~isempty(map)
        I = ind2rgb(I, map);            % indexed (palette) image -> RGB
    end
    if size(I,3) == 4
        I = I(:,:,1:3);                 % drop alpha / 4th channel
    end
    I = im2uint8(I);                    % logical/uint16/double -> 0..255
end

function file = selectImageFile(defaultFile, dialogTitle)
% Ask the user to pick an image when the sample image is not available.
    imageTypes = {'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', ...
                  'Image files (*.png, *.jpg, *.tif, *.bmp)'; ...
                  '*.*', 'All files (*.*)'};
    [fileName, folder] = uigetfile(imageTypes, ...
                                   sprintf('%s (sample %s was not found)', ...
                                           dialogTitle, defaultFile));
    if isequal(fileName, 0)
        error('svd_image_compression:ImageNotSelected', ...
              ['The sample image "%s" is not available. Select an image ', ...
               'when prompted, or place that file in the images/ folder.'], ...
              defaultFile);
    end
    file = fullfile(folder, fileName);
end
