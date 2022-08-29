function [] = SHOW2DMESH(Nds,Tri_Els,Quad_Els,Cface,Cid,ts, varargin)
%SHOW2DMESH Depicts the 2D mesh composed of Tri and Quad Elements
% USAGE:
%	[] = SHOW2DMESH(Nds,Tri_Els,Quad_Els,Cface,Cid,ts);
% INPUTS:
%   Nds		: Nnx2 node locations
%   Tri_Els	: Ntx4 ordered list of nodes in Tri Elements
%		  [Elid n1 n2 n3]
%   Quad_Els	: Nqx5 ordered list of nodes in Quad Elements
%		  [Elid n1 n2 n3 n4]    
%   Cface	: 1x3 Element face color
%            (or) Nelx(1 or 3) Element face values/colors
%   Cid		: 1x3 ID color (negative to switch off)
%   ts		: ID offset (switched off if exactly -100)
% OUTPUTS:
%   Nil
    
    chkfinflag = 0;
    if length(varargin)==1
        chkfinflag=varargin{1};
    end

    Nn = size(Nds,1);
    Nt = size(Tri_Els,1);
    Nq = size(Quad_Els,1);
    Ne = Nt+Nq;
    
    % Plot Tri elements
    for i=1:Nt
        ei = Tri_Els(i,1);
        V = Nds(Tri_Els(i,2:end),:);
        if size(Cface,1)==Nn
            fill(V(:,1),V(:,2),Cface(Tri_Els(i,2:end),:));
        elseif size(Cface,1)==Ne
            if (isfinite(sum(log10(Cface(Tri_Els(i,1),:)))) || chkfinflag==0)
                fill(V(:,1),V(:,2),Cface(Tri_Els(i,1),:))
            else
                plot(V([1:4 1],1), V([1:4 1],2), 'k-')
            end
        else
            fill(V(:,1),V(:,2),Cface(:,1));
        end
        hold on
        
        if isempty(find(Cid<0, 1))
            text(mean(V(:,1))-1e-3*(i/Nt)^(1/3),mean(V(:,2)), ...
                 num2str(ei),'Color',Cid);
        end
    end
    % Plot Quad elements
    for i=1:Nq
        ei = Quad_Els(i,1);
        V = Nds(Quad_Els(i,2:end),:);
        if size(Cface,1)==Nn
            fill(V(:,1),V(:,2),Cface(Quad_Els(i,2:end),:));
        elseif size(Cface,1)==Ne
            if (isfinite(sum(log10(Cface(Quad_Els(i,1),:)))) || chkfinflag==0)
                if(size(Cface, 2) == 1)
                    patch(V(:,1),V(:,2),Cface(Quad_Els(i,1),:))
                else
                    fill(V(:,1),V(:,2),Cface(Quad_Els(i,1),:))
                end
            else
                plot(V([1:4 1],1), V([1:4 1],2), 'k-')
            end
        else
            fill(V(:,1),V(:,2),Cface(1,:));
        end
        hold on
        
        if isempty(find(Cid<0))        
            text(mean(V(:,1))-1e-3*(i/Nq)^(1/3),mean(V(:,2)), ...
                 num2str(ei),'Color',Cid);
        end
    end
    % Plot nodes
    lbls = cellstr(num2str((1:Nn)'));
    if ts~=-100
        text(Nds(:,1),Nds(:,2)+ts*0.2e-3,lbls);%,'Color',[0 0 0]);
        plot(Nds(:,1),Nds(:,2),'ko','MarkerFaceColor','w');
    end
    xlabel('X')
    ylabel('Y')
end