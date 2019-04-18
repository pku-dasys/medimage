function [psnr_Value, mse]=PSNR(A,B)
    maxValue = double(max(A(:)));

    % Calculate MSE, mean square error.
    mseImage = (double(A) - double(B)) .^ 2;
    [rows columns] = size(A);
    
    mse = sum(mseImage(:)) / (rows * columns);

    % Calculate PSNR (Peak Signal to noise ratio)
    psnr_Value = 10 * log10( 256^2 / mse);

end