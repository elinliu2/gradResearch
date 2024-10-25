function tens = conts(tensors,modes,perm)
% only possible to first multiply by an arbitrary number of
% 2-tensors and then multiply by one 3-tensor (other orders will fail!)    
    
% tensors and modes are structs
% tensor contains the tensors to be contracted (multiplied) in order
% modes gives the modes of contraction
% modes{i} = [mode_of_initial_tensor mode_of_next_tensor]
% perm gives the permutation to apply to the final result

% for debugging:
    % flag = 0;
    % if length(tensors) == 4
    %     sA = size(tensors{1});
    %     sB = size(tensors{2});
    %     sC = size(tensors{3});
    %     sD = size(tensors{4});
    %     if length(sA) == 4 && length(sB) == 2 && length(sC) == 2 ...
    %             && length(sD) == 2
    %         if all(sA == [1 1 1 3]) && all(sB == [1 4]) && ...
    %             all(sC == [1 4]) && all(sD == [3 4])
    %             flag = 1;
    %             disp('flag is 1');
    %             sA
    %             sB
    %             sC
    %             sD                
    %         end
    %     end
    % end

   tens = tensors{1};
   num_modes = length(size(tens));
   current_shift = 0;
   
   for i=1:length(modes)
       next = tensors{i+1};

       mode1 = modes{i}(2);
       
       mode2 = mod(modes{i}(1)-current_shift,num_modes);
       if mode2 == 0
           mode2 = num_modes;
       end

       % for debugging:
       % if flag
       %    i
       %    size(next)
       %    size(tens)
       %    mode1
       %    mode2
       % end
       
       tens = con2(next,tens,mode1,mode2,[]);
       %tens = squeeze(tens);
       
       current_shift = current_shift + mode2-1;
       current_shift = mod(current_shift,num_modes);

       % for debugging:
       % if flag
       %    current_shift
       %    size(tens)
       % end
   end

   if length(size(next)) == 2
       tens = shiftdim(tens,mod(-current_shift,num_modes));
   elseif length(size(next)) == 3
       tens = shiftdim(tens,mod(-current_shift,num_modes)+1);
   end
   if ~isempty(perm)   
       tens = permute(tens,perm);
   end
   
   tens = squeeze(tens);
   
   if isempty(tens)
       if size(tens,1) == 0
           return;
       elseif size(tens,2) == 0
           tens = shiftdim(tens,1);
       elseif size(tens,3) == 0
           tens = shiftdim(tens,2);
       elseif size(tens,4) == 0
           tens = shiftdim(tens,3);
       end
   end

end
