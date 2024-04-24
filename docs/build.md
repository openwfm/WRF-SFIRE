# For Contributors to the Docs

## Build and view the docs locally 

### Install mkdocs

     conda create -n mkdocs -c conda-forge mkdocs mkdocstrings 
     conda activate mkdocs
     pip install markdown_include

### Clone this repo and build the docs

    git clone https://github.com/openwfm/WRF-SFIRE
    cd WRF-SFIRE
    conda activate mkdocs
    mkdocs build

** All files paths below are relative to the directory WRF-SFIRE **

### View the site

Open the file `docs/index.html` in a web browser.

##  How to write the docs

### Make your own branch

If you are contributing code, make your own branch as you normaly would per the instructions elsewhere
(usually off the develop or release branch) and make changes in the documentation
to document your changed code.

When you are contributing to the docs only and do not change the code, make a new branch off the `master`
branch and give it a name indicating this is a documentation change, such as `docs-`something. 

### Edit the docs

Chapters are Markdown files in directory `docs`. Please stay close to the style of the files already there.
The navigation bar on the left (or a pancake when viewing the site on a small screen) is created from the
`nav:` section at the bottom of file `mkdocs.yml`. When you add a new Markdown file, add also a corresponding line there.


### Graphics

You can link to graphics staged elsewhere. If you *must* have graphics files, you can put *small* graphics files only 
in `docs/resources` and link to them locally.  Reduce the resolution of images and do not commit any high-resolution large files.
Binary files like images  bloat the repository  size forever, even of they are replaced or deleted in later commits.

### Submit your changes

Merge the branch you started from into your branch and resolve any conflicts. Then test the result. Merge again in case anything 
changed in the meantime. If nothing changed, you can make a Github pull request to merge your changes in.






## Resources

https://www.markdownguide.org



[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

