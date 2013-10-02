classdef TypeHierarchy < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tags=cell(0,0);
        rmap=zeros(0,0);
    end
    
    methods
        function [obj] = TypeHierarchy()
            obj.tags={'root'};
            obj.rmap=[1,0];
            return;
        end
        
        function [ok] = addAll(obj,parent, node_names)
            ok=0;
            if ischar(parent)==0
                return;
            end
            
            for i=1:length(node_names)
                obj.addNode(parent, node_names{i});
            end
            ok=1;
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
                node_name_id = length(obj.tags);

                obj.rmap = [obj.rmap; parent_id, node_name_id];
                ok=1;
                return;
            end
        end
        
        function [ok]=removeNode(obj,node_name)
            ok=0;
            if ischar(node_name)==0 || strcmp(node_name,'root')==1
                return;
            end
            
            if sum(ismember(obj.tags, node_name))==1
                node_id = find(ismember(obj.tags, node_name)==1)
                obj.rmap(find(obj.rmap(:,1)==node_id | obj.rmap(:,2)==node_id),:)=[];
                
                %don't remove the tag
                %obj.tags{node_id,1}='';
            end     
        end
        
        function [parent] = getParent(obj,node_name)
            parent =-1;
            if ischar(node_name)==0
                return;
            end
            
            if sum(ismember(obj.tags, node_name))==1
                node_name_id = find(ismember(obj.tags, node_name)==1);
                id = unique(obj.rmap(find(obj.rmap(:,2)==node_name_id),1));
                if ~isempty(id)
                    parent = obj.tags{id(1,1),1};
                end
            end
            
        end
        
        function children = getChildren(parent)
            children = -1;
            if ischar(parent)==0
                return;
            end
            
            parent_id = find(ismember(obj.tags, parent)==1);
            children_ids = obj.rmap(find(obj.ramp(:,1)==parent_id),2);
            children = obj.tags{children_ids};
        end
        
        function [ok]=printHierarchy()
        end
    end
    
end

