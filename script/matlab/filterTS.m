% Y(vertex,time)

function [Y_filtered] = filterTS(Y, t, interval)
    [nrow,ncol] = size(Y);
    Y_filtered = zeros(nrow, ncol);
    for i = 1:nrow
        y_TS = timeseries(Y(i,:));
        y_filtered_TS = idealfilter(y_TS,interval,'pass');
        Y_filtered(i,:) = reshape(y_filtered_TS.Data,1,size(Y,2));
    end
end