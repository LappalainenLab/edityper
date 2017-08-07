#!/usr/bin/env python2

'''Configure the CRISPR program'''

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this moduel: " + __name__)


import os
import time
import logging
import ConfigParser

try:
    from write_sam import READ_GROUP_VALUES
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)


SECTION_NAME = 'CRISPR' # type: str
READ_SECTION = 'READ_GROUP' # type: str
CONFIG_VALUES = { # type: Dict[str, str]
    # Option name       : option type
    'analysis_mode'     :   'str',
    'outdirectory'      :   'str',
    'project'           :   'str',
    'input_file'        :   'str',
    'sample_list'       :   'str',
    'reference'         :   'str',
    'template'          :   'str',
    'gap_opening'       :   'int',
    'gap_extension'     :   'int',
    'pvalue_threshold'  :   'float'
}

def _validate_config(conf_dict, typecheck=False): # type: (dict, bool) -> None
    logging.info("Validating configuration data")
    validation_start = time.time() # type: float
    required = { # type: Set[str]
        'analysis_mode',
        'gap_opening',
        'gap_extension',
        'pvalue_threshold',
        'reference',
        'template',
        'outdirectory',
        'project'
    }
    #   Check the required values
    for value in required: # type: str
        if value not in conf_dict:
            raise ValueError("Could not find option '%s' in config file" % value)
        if typecheck and not isinstance(conf_dict[value], eval(CONFIG_VALUES[value])):
            raise ValueError("Option %s must be of type %s" % (value, CONFIG_VALUES[value]))
    try:
        assert ('input_file' in conf_dict and 'sample_list' not in conf_dict) or \
            ('sample_list' in conf_dict and 'input_file' not in conf_dict)
    except AssertionError:
        raise ValueError("Cannot have both an input file and sample list")
    logging.debug("Configuration validation took %s seconds", round(time.time() - validation_start, 3))


def make_prefix_name(directory, base): # type: (str, str) -> str
    """Combine a directory and base together"""
    directory = os.path.abspath(directory)
    if directory[-1] != '/':
        directory += '/'
    return ''.join((directory, base))


def mkdir(directory): # type: (str) -> None
    """Ensure a directory exists. If not, make it and any parents that do not yet exist"""
    if not os.path.isdir(directory):
        os.makedirs(directory)


def write_config(args): # type: (Dict[str, Any]) -> None
    """Write a configuration file"""
    _validate_config(conf_dict=args)
    try:
        cfile_name = make_prefix_name(directory=args['outdirectory'], base=args['project']) + '.ini'
        mkdir(directory=args['outdirectory'])
        logging.info("Creating config file")
        write_start = time.time() # type: float
    except KeyError:
        sys.exit(logging.critical("We need a config file to write to"))
    config = ConfigParser.ConfigParser() # type: ConfigParser.ConfigParser
    config.add_section(SECTION_NAME)
    for option in sorted(CONFIG_VALUES): # type: str
        try:
            config.set(SECTION_NAME, option, str(args[option]))
        except KeyError:
            logging.warning("Option '%s' not found in argument dictionary, skipping...", option)
            continue
        else:
            logging.info("Setting option '%s' with value '%s'", option, args[option])
    #   read group configuration
    if filter(None, ('read' in arg for arg in args)):
        logging.info("Setting read group options")
        config.add_section(READ_SECTION)
        for option in sorted(READ_GROUP_VALUES):
            try:
                config.set(READ_SECTION, option, str(args[option]))
            except KeyError:
                logging.warning("Read group option '%s' not found in argument dictionary, skipping...'", option)
                continue
            else:
                logging.info("Setting option '%s' with value '%s'", option, args[option])
    with open(cfile_name, 'w') as config_file:
        config.write(config_file)
        logging.info("Writing config to config file %s", cfile_name)
    logging.debug("Writing configuration file took %s seconds", round(time.time() - write_start, 3))


def read_config(config_file): # type: (str) -> Dict[str, Any]
    """Read a configuration file"""
    logging.info("Using configuration file %s", config_file)
    read_start = time.time() # type: float
    config = ConfigParser.ConfigParser() # type: ConfigParser.ConfigParser
    config.readfp(open(config_file, 'r'))
    conf_d = dict() # type: Dict[str, Any]
    #   Load the values from the config file as their specific type
    for option, opt_type in CONFIG_VALUES.items(): # type: str, str
        try:
            if opt_type == 'int':
                conf_d[option] = config.getint(SECTION_NAME, option)
            elif opt_type == 'float':
                conf_d[option] = config.getfloat(SECTION_NAME, option)
            elif opt_type == 'bool':
                conf_d[option] = config.getboolean(SECTION_NAME, option)
            else:
                conf_d[option] = config.get(SECTION_NAME, option)
        except ConfigParser.NoOptionError: # If we don't have this option, we don't particularly care
            logging.warning("Tried setting option '%s', but option was not found in config file", option)
            continue
        else:
            logging.info("Option '%s' set to '%s'", option, conf_d[option])
    for option in READ_GROUP_VALUES:
        try:
            conf_d[option] = config.get(READ_SECTION, option)
        except ConfigParser.NoOptionError: # If we don't have this option, we don't particularly care
            logging.warning("Tried setting read group option '%s', but option was not found in config file", option)
            continue
        except ConfigParser.NoSectionError:
            logging.warning("No read group options found in config file, will make minimal read group headers")
            break
        else:
            logging.info("Read group option '%s' set to '%s'", option, conf_d[option])
    logging.debug("Loading the config file took %s seconds", round(time.time() - read_start, 3))
    _validate_config(conf_dict=conf_d, typecheck=True)
    mkdir(conf_d['outdirectory'])
    return conf_d
