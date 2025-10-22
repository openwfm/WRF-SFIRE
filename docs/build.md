# For Editors of the Docs

## Files

    git clone https://github.com/openwfm/WRF-SFIRE
    cd WRF-SFIRE

    .readthedocs.yaml`    # ReadTheDocs configuration file
    mkdocs.yml            # MkDocs configuration file.
    docs/
         requirements.in  # PIP configuration file.
         requirements.txt # PIP configuration file.
         index.md         # The documentation homepage.
         ...              # Other markdown pages, images and other files.

## Build the software 

The docs in the `master` branch of 
[https://github.com/openwfm/WRF-SFIRE](https://github.com/openwfm/WRF-SFIRE)
are displayed at [https://wrf-sfire.readthedocs.io](https://wrf-sfire.readthedocs.io) 

You can preview the docs in your own copy or in other branches locally as follows. First install mkdocs (once):

     conda create -n mkdocs python=3.12
     conda activate mkdocs
     pip install -r docs/requirements.txt

## Commands

Every time:

   conda activate mkdocs
   mkdocs build  #  Build the documentation site.

and open the file site/index.html. This is far from perfect, though; if clicking on a link gives you a directory,
open file index.html in that directory.

Other commands available:

  mkdocs serve # Start the live-reloading docs server.
  mkdocs -h    #  Print help message and exit.

For full documentation visit [mkdocs.org](https://www.mkdocs.org).

##  How to edit 

### Make your own branch

If you are contributing code, make your own branch as you normaly would per the instructions elsewhere
(usually off the develop or release branch) and make changes in the documentation
to document your changed code.

When you are contributing to the docs only and do not change the code, make a new branch off the `master`
branch and give it a name indicating this is a documentation change, such as `docs-`something. 

### Edit the docs

Chapters are Markdown files in directory `docs`.
The navigation bar on the left (or a pancake when viewing the site on a small screen) is created from the
`nav:` section at the bottom of file `mkdocs.yml`. When you add a new Markdown file, you need to add also a corresponding line there.

### Graphics

Please reduce the resolution of images and do not commit any large files.
Large binary files like images bloat the repository size forever, 
even of they are replaced or deleted in later commits.

### Submit your changes

Merge the branch you started from into your branch and resolve any conflicts. 
Then test the result. Merge again in case anything 
changed in the meantime. If nothing changed, you can make a pull request 
on Github.

## Resources

https://www.markdownguide.org

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

