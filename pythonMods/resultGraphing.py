__author__ = 'sam'


import json
import outputProcessing
import matplotlib
import numpy as np
import mhdTools
import os
from itertools import imap, izip
import graph_tool.all as gt


class StructureStatisticsMap(outputProcessing.StructureFinder):

    def __init__(self, structure_data_path, metadata_path):
        outputProcessing.StructureFinder.__init__(self, structure_data_path)
        self.metadata = self.load_metadata(metadata_path)
        self.metadataDirPath = os.path.dirname(metadata_path)
        self.exclusionMap = self._initialize_map()
        self.occurrenceMap = self.exclusionMap.copy()
        self.pIncludeMap = self.exclusionMap.copy()
        self._populate_maps()

    def _initialize_map(self):
        return {id_no: 0 for id_no in self._generate_structure_ids(self.structureData)}

    def _populate_maps(self):
        # Load occurrence data from images
        self._load_occurrence_data()
        self._load_exclusion_data()
        self._compute_inclusion_probability()

    def _compute_inclusion_probability(self):
        self._accumulate_inclusion_probability(self.structureData)

    def _accumulate_inclusion_probability(self, node):
        node_id = node['id']
        if node['children']:
            map(self._accumulate_inclusion_probability, node['children'])
            self.pIncludeMap[node_id] = np.mean(np.array([self.pIncludeMap[child['id']] for child in node['children']]))
        elif self.occurrenceMap[node_id] > 0:
            p_exclude = float(self.exclusionMap[node_id])/float(self.occurrenceMap[node_id])
            self.pIncludeMap[node_id] = 1.0 - p_exclude
        else:
            self.pIncludeMap[node_id] = 0.0

    def _load_exclusion_data(self):
        for image_data in self.metadata:
            # TODO: There are 3 separate things that the matlab JSON api will put into this list.
            #  None if no exclusions
            #  List of ints if one exclusion
            #  List of lists if there are > 1 exclusion
            # This is ugly and must be fixed
            if image_data['regionIdsToExclude'] is None:
                pass
            elif type(image_data['regionIdsToExclude'][0]) is int:
                id_no = image_data['regionIdsToExclude'][0]
                if id_no > 0:
                    self.exclusionMap[id_no] += 1
            elif type(image_data['regionIdsToExclude'][0]) is list:
                for id_no, _ in image_data['regionIdsToExclude']:
                    # Exclusions are specific to a hemisphere, Ids will be listed twice if they are excluded from both
                    if id_no > 0:
                        self.exclusionMap[id_no] += 1
            else:
                pass

    def _load_occurrence_data(self):
        for image_data in self.metadata:
            if not image_data['sliceUsable'] == 'no':
                annotation_image_name = os.path.basename(image_data['registeredAtlasLabelsPath'])
                annotation_image_path = os.path.join(self.metadataDirPath, annotation_image_name)
                annotation_image, _ = mhdTools.load_mhd(annotation_image_path)

                included_ids = np.unique(annotation_image)
                # Background is zero, we don't care about the background
                included_ids = list(included_ids[included_ids > 0])
                self._accumulate_ids(included_ids, self.structureData)
                for id_no in included_ids:
                    # If a structure is present in the labeled image, it is present exactly twice
                    self.occurrenceMap[id_no] += 2

    def _accumulate_ids(self, included_ids, node):
        for child in node['children']:
            self._accumulate_ids(included_ids, child)
        if any((child['id'] in reversed(included_ids) for child in node['children'])):
            included_ids.append(node['id'])

    @staticmethod
    def load_metadata(metadata_path):
        metadata_file = open(metadata_path)
        metadata = json.load(metadata_file)
        metadata_file.close()
        return metadata


class StructureStatisticsGrapher():

    def __init__(self, root_node=None):
        self.structureGraph = gt.Graph(directed=True)
        self.vertexIds = self.structureGraph.new_vertex_property('int32_t')
        if root_node is not None:
            self.add_with_children(root_node)

    def add_with_children(self, node):
        """
        Adds the given structure node and all of it's children to the structureGraph.
        Each child node is linked to it's parent

        The structure node corresponding to each new vertex is saved in the vertexIds propertyMap

        :param node: 'root' node to be added to the graph
        :return: root node's vertex
        """
        vertex = self.structureGraph.add_vertex()
        self.vertexIds[vertex] = node['id']
        for child_vertex in imap(self.add_with_children, node['children']):
            self.structureGraph.add_edge(child_vertex, vertex)

        return vertex

    def add_property_with_id_map(self, property_map, property_type):
        prop = self.structureGraph.new_vertex_property(property_type)
        for vert in self.structureGraph.vertices():
            id_no = self.vertexIds[vert]
            if id_no in property_map.keys():
                prop[vert] = property_map[id_no]

        return prop

    def draw_radial_graph(self,
                          output_path,
                          vertex_text=None,
                          vertex_fill_color=[0.640625, 0, 0, 0.9],
                          bg_color=(.75, .75, .75, 1),
                          output_size=(1000, 1000),
                          vcmap=matplotlib.cm.summer,
                          **kwargs):
        """
        Draws the structureGraph as a radial graph to output_path
        :param output_path: path to save the image of the graph to
        :return:
        """
        pos = gt.radial_tree_layout(self.structureGraph, self.structureGraph.vertex(0))
        gt.graph_draw(self.structureGraph,
                      pos=pos,
                      output=output_path,
                      bg_color=bg_color,
                      vertex_text=vertex_text,
                      output_size=output_size,
                      vertex_fill_color=vertex_fill_color,
                      vcmap=vcmap,
                      **kwargs)


def test_stats_map():
    metadata_path = '/home/sam/Desktop/imagesForJin_modifiedAtlas/metadata.json'
    structure_data_path = '/home/sam/code/fishRegistration/structureData.json'
    stats_map = StructureStatisticsMap(structure_data_path, metadata_path)
    acr_gen_fcn = stats_map.get_structure_property_generator_function('acronym')
    id_gen_fcn = stats_map.get_structure_property_generator_function('id')
    for acr, id_no in izip(acr_gen_fcn(stats_map.structureData), id_gen_fcn(stats_map.structureData)):
        print acr + ': {} {} {}'.format(stats_map.occurrenceMap[id_no], stats_map.exclusionMap[id_no], stats_map.pIncludeMap[id_no])
    return stats_map


def test_structure_graph():
    stats_map = test_stats_map()
    structure_data_path = '/home/sam/code/fishRegistration/structureData.json'
    grapher = StructureStatisticsGrapher(stats_map.structureData)

    acr_gen_fcn = stats_map.get_structure_property_generator_function('acronym')
    id_gen_fcn = stats_map.get_structure_property_generator_function('id')
    id_list = list(id_gen_fcn(stats_map.structureData))
    acronym_map = {id_no: acr for acr, id_no in izip(acr_gen_fcn(stats_map.structureData), id_list)}
    acronym_property = grapher.add_property_with_id_map(acronym_map, 'string')

    colormap = generate_colormap(stats_map.pIncludeMap)
    color_property = grapher.add_property_with_id_map(colormap, 'double')

    grapher.draw_radial_graph('/home/sam/Desktop/test_graph.png',
                              vertex_text=acronym_property,
                              vertex_fill_color=color_property)


def get_absent_regions_filter_dict(stats_map):
    annotation_atlas, _ = mhdTools.load_mhd('/home/sam/code/fishRegistration/atlasVolume/annotation.mhd')
    annotation_inclusion_list = list(np.unique(annotation_atlas))
    stats_map._accumulate_ids(annotation_inclusion_list, stats_map.structureData)

    filter_dict = dict()
    id_gen_fcn = stats_map.get_structure_property_generator_function('id')
    for id_no in id_gen_fcn(stats_map.structureData):
        # id is included if True
        filter_dict[id_no] = id_no in annotation_inclusion_list

    print all(list(filter_dict.values()))
    print any(list(filter_dict.values()))

    return filter_dict


def get_filtered_structure_data(filter_map, structure_data):
    filtered_data = structure_data.copy()
    filtered_data['children'] = get_filtered_children(filter_map, structure_data['children'])
    return filtered_data


def get_filtered_children(filter_map, children):
    filtered_children = list()
    for child in children:
        if filter_map[child['id']]:
            filtered_children.append(get_filtered_structure_data(filter_map, child))
        else:
            print child['name']

    return filtered_children


def create_figures():
    metadata_path = '/home/sam/Desktop/imagesForJin_modifiedAtlas/metadata.json'
    # structure_data_path = '/home/sam/code/fishRegistration/structureData.json'
    structure_data_path = '/home/sam/code/fishRegistration/structureData_filtered.json'
    stats_map = StructureStatisticsMap(structure_data_path, metadata_path)

    # The stats map is also a structure finder, this is probably bad practice but I'm feeling very lazy...
    acr_gen_fcn = stats_map.get_structure_property_generator_function('acronym')
    id_gen_fcn = stats_map.get_structure_property_generator_function('id')

    # Some different nodes that we are interested in
    root_node = stats_map.structureData  # full data set
    isocortex_node = stats_map.search_structure_data_for_attribute('name', 'Isocortex')
    AUD_node = stats_map.search_structure_data_for_attribute('acronym', 'AUD')
    TEa_node = stats_map.search_structure_data_for_attribute('acronym', 'TEa')
    ECT_node = stats_map.search_structure_data_for_attribute('acronym', 'ECT')
    PERI_node = stats_map.search_structure_data_for_attribute('acronym', 'PERI')

    # Cortical Subplate contains amygdala, Claustrum and LA
    CTXsp_node = stats_map.search_structure_data_for_attribute('acronym', 'CTXsp')

    # Filter dict to filter nodes that don't appear in the brain atlas
    filter_map = get_absent_regions_filter_dict(stats_map)
    filtered_node = get_filtered_structure_data(filter_map, stats_map.structureData)

    # Colormap based on 1 - pInclude
    colormap = generate_colormap(stats_map.pIncludeMap)
    edgecolor = generate_edge_cm(stats_map.occurrenceMap)

    # mapping of id_no : acronyms for every structure
    acronym_map = {id_no: acr for acr, id_no in izip(acr_gen_fcn(root_node), id_gen_fcn(root_node))}


    # Graph for the full data set
    full_graph = StructureStatisticsGrapher(root_node)
    ecm_property = full_graph.add_property_with_id_map(edgecolor, 'double')
    cm_property = full_graph.add_property_with_id_map(colormap, 'double')
    acr_property = full_graph.add_property_with_id_map(acronym_map, 'string')
    full_graph.draw_radial_graph('/home/sam/Desktop/full_graph.png',
                                 output_size=(4500, 4500),
                                 vertex_text=acr_property,
                                 vertex_color=ecm_property,
                                 vertex_fill_color=cm_property)


    ## Filtered data set
    #filtered_graph = StructureStatisticsGrapher(filtered_node)
    #ecm_property = filtered_graph.add_property_with_id_map(edgecolor, 'double')
    #cm_property = filtered_graph.add_property_with_id_map(colormap, 'double')
    #acr_property = filtered_graph.add_property_with_id_map(acronym_map, 'string')
    #filtered_graph.draw_radial_graph('/home/sam/Desktop/filtered_graph.png',
    #                                 output_size=(4500, 4500),
    #                                 vertex_text=acr_property,
    #                                 vertex_color=ecm_property,
    #                                 vertex_fill_color=cm_property)

    # Graph for the Isocortex
    iso_graph = StructureStatisticsGrapher(isocortex_node)
    ecm_property = iso_graph.add_property_with_id_map(edgecolor, 'double')
    cm_property = iso_graph.add_property_with_id_map(colormap, 'double')
    acr_property = iso_graph.add_property_with_id_map(acronym_map, 'string')
    iso_graph.draw_radial_graph('/home/sam/Desktop/iso_graph.png',
                                output_size=(4000, 4000),
                                vertex_text=acr_property,
                                vertex_color=ecm_property,
                                vertex_fill_color=cm_property)

    # Graph for the Auditory Cortex
    ACtx_graph = StructureStatisticsGrapher()
    ecm_property = ACtx_graph.add_property_with_id_map(edgecolor, 'double')
    v0 = ACtx_graph.structureGraph.add_vertex()
    for node in (AUD_node, TEa_node, ECT_node, PERI_node):
        new_vertex = ACtx_graph.add_with_children(node)
        ACtx_graph.structureGraph.add_edge(new_vertex, v0)
    cm_property = ACtx_graph.add_property_with_id_map(colormap, 'double')
    acr_property = ACtx_graph.add_property_with_id_map(acronym_map, 'string')
    ACtx_graph.draw_radial_graph('/home/sam/Desktop/auditory_cortex_graph.png',
                                 output_size=(2000, 2000),
                                 vertex_font_size=24,
                                 vertex_text=acr_property,
                                 vertex_color=ecm_property,
                                 vertex_fill_color=cm_property,
                                 vertex_size=50)

    # Graph for the Cortical Subplate
    ctxsp_graph = StructureStatisticsGrapher(CTXsp_node)
    ecm_property = ctxsp_graph.add_property_with_id_map(edgecolor, 'double')
    cm_property = ctxsp_graph.add_property_with_id_map(colormap, 'double')
    acr_property = ctxsp_graph.add_property_with_id_map(acronym_map, 'string')
    ctxsp_graph.draw_radial_graph('/home/sam/Desktop/cortical_subplate_graph.png',
                                  output_size=(1750, 1750),
                                  vertex_text=acr_property,
                                  vertex_color=ecm_property,
                                  vertex_fill_color=cm_property,
                                  vertex_font_size=36)


def generate_colormap(p_map):
    colormap = dict()
    for id_no, p in p_map.iteritems():
        colormap[id_no] = 1 - p

    return colormap


def generate_edge_cm(o_map):
    return {id_no: num_occur < 5 for id_no, num_occur in o_map.iteritems()}


def main():
    create_figures()

if __name__ == '__main__':
    main()
