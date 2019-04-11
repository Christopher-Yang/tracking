function [e, a] = error_unwrap(err, avg, uw)

    switch uw
        case {1,2,4,5,9,10,11,12,14,15,16,17,18,19,20,22,23,24}
%         case {1, 2, 4, 5}
            e = err.*(180/pi);
            a = avg;
        case 3
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
        case 6
            e = err;
%             e(2,7,end-25) = e(2,7,end-25) - 2*pi;
            
            e = e.*(180/pi);
            
            a = avg;
            a(2,7) = a(2,7) - 2*pi;
        case 7
            e = err;
%             e(2,1,26) = e(2,1,26) - pi;
%             e(3,1,26) = e(3,1,26) + pi;
%             e(3,2:7,26) = e(3,2:7,26) + 2*pi;
%             
%             e(2:3,:,end-25) = e(2:3,:,end-25) - pi;
%             e(4:5,2:7,end-25) = e(4:5,2:7,end-25) - pi;
            e = e.*(180/pi);
            
            a = avg;
            a(3,:) = a(3,:) + 2*pi;
%             a(2:3,1) = a(2:3,1) - pi;
        case 8
            e = err;
% %             e(:,:,end-25) = e(:,:,end-25) - pi;
            e = e.*(180/pi);
%             
            a = avg;
        case 13
            e = err;
            e(2,6,end-25) = e(2,6,end-25) - 2*pi;
            e = e.*(180/pi);
            
            a = avg;
            a(2,5:7) = a(2,5:7) - 2*pi;
        case 21
            e = err;
            e(2,6:7,end-25) = e(2,6:7,end-25) + 2*pi;
            e(2,6:7,26) = e(2,6:7,26) + 2*pi;
            e = e.*(180/pi);
            
            a = avg;
    end
end