function [R] = findzeroQ(qa)

   % return a quaternion matrix with entries 0/1 indicating the zero
   % entries of qa

    w = (qa.w==0); 
    x = (qa.x == 0);
    y = (qa.y == 0);
    z = (qa.z == 0);

    R = quaternion(w,x,y,z);

end

