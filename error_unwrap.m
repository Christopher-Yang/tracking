function [e, a] = error_unwrap(err, avg, uw)

    switch uw
        case {1,2,3,4,5,6,7,8,9,10,11,12}
            
            e = err;
            a = avg;
            
            e = e.*(180/pi);
        case 13
            e = err;
            e(1:5,:,26) = e(1:5,:,26) + pi;
            e(1,2,end-25) = e(1,2,end-25) - pi;
            e(3,1,end-25) = e(3,1,end-25) - pi;
            e(3,2:7,end-25) = e(3,2:7,end-25) + pi;
            e([2 4:5],:,end-25) = e([2 4:5],:,end-25) + pi;
            e = e.*(180/pi);
            
            a = avg;
            a(2,:) = a(2,:) - pi;
            a(3:5,:) = a(3:5,:) + pi;
    end
end