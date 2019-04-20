function spectrogram = spectrogram_plus(x, fs, fft_length, window_length, overlap)
    
    
    sx = size(x)
    if (sx(1) == 1)
        x = x';
    end

    numFFTs = floor(length(x)/(window_length-overlap));
    
    fft_mat = zeros(fft_length, numFFTs);
    
    win = hann(window_length);
    
    for ii = 1:numFFTs
        
        fft_mat(:, ii) = fft_plus(win.*x((ii-1)*(window_length-overlap)+1:(ii-1)*(window_length-overlap)+window_length), fs, fft_length);
        
    end

    spectrogram = fft_mat((fft_length/2+1):end,:);

end