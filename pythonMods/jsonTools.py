import json
import os


class MetadataHandler():
    """
    A class to handle the metadata.json file used in the image registration pipeling

    I handle loading, and updating the data in metadata.json
    """

    metadata = None

    def __init__(self, experiment_path=os.getcwd()):
        self.experimentPath = experiment_path
        self.metadataPath = self.generate_metadata_path(experiment_path)

    def load_metadata(self):
        """
        Loads the metadata from metadata.json and stores it at self.json
        :return: Pointer to self.metadata
        """
        json_file = open(self.metadataPath, 'r')
        self.metadata = json.load(json_file)
        json_file.close()
        return self.metadata

    def update_metadata(self):
        """
        Updates the data stored at metadata.json with the data in self.metadata
        :return: zilch
        """
        json_file = open(self.metadataPath, 'w')
        json.dump(self.metadata, json_file, sort_keys=True, indent=4)
        json_file.close()

    def ensure_metadata_directory(self):
        """
        Creates self.experimentPath/.registrationData if it doesn't exist
        :return: nope
        """
        dir_path = os.path.join(self.experimentPath, '.registrationData')
        _ensure_dir(dir_path)

    @staticmethod
    def generate_metadata_path(experiment_path):
        """
        Generates the path the the metadata.json file
        :param experiment_path: root path for the experiment data.  Must contain the .registrationData directory
        :return: a string containing the path to the .metadata.json file: 'experimentPath/.registrationData.metadata.json'
        """
        return os.path.join(experiment_path, '.registrationData', 'metadata.json')


def _ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    elif not os.path.isdir(dir_path):
        print "Warning! Non-directory file exists at:\n\t%s" % dir_path

