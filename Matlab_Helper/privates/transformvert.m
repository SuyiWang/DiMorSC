function vert = transformvert(vert, trans, scale)
for i=1:3
    vert(:,i) = vert(:,i).*scale(i) + trans(i);
end