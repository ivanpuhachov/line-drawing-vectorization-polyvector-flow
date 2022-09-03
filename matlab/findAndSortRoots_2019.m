function [ curRoots ] = findAndSortRoots_2019( Xi )
curRoots = roots([1 0 Xi(2) 0 Xi(1)]);
[~,I] = sort(abs(curRoots),'descend');
curRoots = curRoots(I);
[~,I] = sort(angle(curRoots(1:2)));
curRoots(1:2) = curRoots(I);
[~,I] = sort(angle(curRoots(3:4)));
curRoots(3:4) = curRoots(I+2);
end

