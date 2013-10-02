classdef TypeHierarchy < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tags=cell{0,0};
        rmap=zeros(0,0);
    end
    
    methods
        function [obj] = TypeHierarchy()
        end
        
        function [ok]=addNode(obj,parent, node_name)
            ok=0;
            
            % both must be strings
            if ischar(parent)==0 || ischar(node_name)==0
                return;
            end
            
            if sum(ismember(obj.tags, parent))==1 && sum(ismember(obj.tags, node_name))==1
                parent_id = find(ismember(obj.tags, parent)==1);
                node_name_id = find(ismember(obj.tags, node_name)==1);
                
                % if it is not already an entry in the rmap
                if isempty(intersect(find(obj.rmap(:,1)==parent_id), find(obj.rmap(:,2)==node_name_id)))
                    obj.rmap = [obj.rmap; parent_id, node_name_id];
                    ok=1;
                    return;
                end
            elseif sum(ismember(obj.tags, parent))==1 && sum(ismember(obj.tags, node_name))==0
                parent_id = find(ismember(obj.tags, parent)==1);
                obj.tags{length(obj.tags)+1,1} = node_name;
                node_name_id = length(obj.tags)+1;

                obj.rmap = [obj.rmap; parent_id, node_name_id];
                ok=1;
                return;
            end
        end
        
        function [ok]=removeNode(node_name)
            ok=0;
            if ischar(node_name)==0
                return;
            end
            
            if sum(ismember(obj.tags, node_name))==1
                node_id = find(ismember(obj.tags, node_name)==1);
                rmap[find(rmap(:,1)==node_id || rmap(:,2)==node_id),:]=[];
                obj.tags(node_id,:)=[];
            end
                
        end
        
        function [parent] = getParent(node_name)
            parent =-1;
            if ischar(node_name)==0
                return;
            end
            
            
            
        end
        
        function [ok]=printHierarchy()
        end
    end
    
end

