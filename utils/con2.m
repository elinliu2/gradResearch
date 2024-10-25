
% contract tensors A and B along Aindex for A and Bindex for B
% using shiftdim; final answer output in order perm


function C = con2(A,B,Aindex,Bindex,perm)

  if Aindex > length(size(A))
     Aindex = length(size(A)); 
  end
  
  if Bindex > length(size(B))
     Bindex = length(size(B)); 
  end    
    
  Alength = size(A,Aindex);
  Blength = size(B,Bindex);
  
  if Alength ~= Blength

      % checking if one or both tensors is a number
      if length(size(A)) == 2
          if all(size(A) == [1 1])
              C = A*B;
              return;
          end
      end      
      if length(size(B)) == 2
          if all(size(B) == [1 1])
             C = A*B;
             return;
          end
      end
      
      %checking if a tensor contraction dimension is one
      if Alength == 1
          B = reshape(B,[1 size(B)]);
          Bindex = 1;
      end
      if Blength == 1
          A = reshape(A,[1 size(A)]);
          Aindex = 1;
      end
      
      Alength = size(A,Aindex);
      Blength = size(B,Bindex);
      
      if Alength ~= Blength      
          disp('Contraction fails due to mismatch in lengths.');
          keyboard
      end
  end
  
  reshaping = [];
  for i=1:length(size(A))-1
      if Aindex+i > length(size(A))
          adding = Aindex+i - length(size(A));
      else
          adding = Aindex+i; 
      end
      reshaping = [reshaping size(A,adding)];
  end
  
  for j=1:length(size(B))-1
      if Bindex+j > length(size(B))
          adding = Bindex+j - length(size(B));
      else
          adding = Bindex+j;
      end
      reshaping = [reshaping size(B,adding)];
  end
  
  if Alength == 0
      tempA = size(shiftdim(A,Aindex));
      tempB = size(shiftdim(B,Bindex-1));
      C = zeros([tempA(1:end-1) tempB(2:end)]);
      return;
  end
        
  if ~isempty(perm)      
      C = permute(reshape(full(reshape(shiftdim(A,Aindex),[],Alength)*...
          sparse(reshape(shiftdim(B,Bindex-1),Blength,[]))),reshaping),perm);
  else
      C = reshape(full(reshape(shiftdim(A,Aindex),[],Alength)*...
          sparse(reshape(shiftdim(B,Bindex-1),Blength,[]))),reshaping);
  end
  
end