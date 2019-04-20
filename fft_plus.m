function [X, f] = fft_plus(x, fs, n)

    % regular fft and normalize
    
    if(n)
        
        X_old = fft(x, n);

        X_old = X_old./n;%'./n;
        size(X_old)
    else
        X_old = fft(x')/length(x); %/length(x);
    end
    
    % shift fft values greater than fs/2
    %X_old
    X = zeros(1, length(X_old));
    %X
    
    %left side
    X(1:ceil(length(X_old)/2)-1) = X_old(floor(length(X_old)/2)+2:end);
    
    %right side
    X(ceil(length(X_old)/2):end) = X_old(1:floor(length(X_old)/2)+1);
    
    % create freq vector
    df = fs/length(x);
    f = df*(-length(X)/2+1:length(X)/2);
    
end