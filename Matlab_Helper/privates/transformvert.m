function vert = transformvert(vert, trans, scale)
for i=1:2
    vert(:,i) = vert(:,i).*scale(i) + trans(i);
end