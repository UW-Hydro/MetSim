Making a new release of Metsim is easy. Just follow the following steps:

Release per project:

*   Update release notes in docs/whats-new.rst

*   Tag commit

        git tag -a x.x.x -m 'Version x.x.x'

*   and push to github

        git push upstream master --tags

*  Upload to PyPI

        git clean -xfd
        python setup.py sdist bdist_wheel --universal
        twine upload dist/*

*   Update conda recipe feedstocks on `conda-forge <https://conda-forge.github.io/>`_.

    *  Update conda-smithy and run conda-smithy rerender

            git clone git@github.com:conda-forge/metsim-feedstock
            cd metsim-feedstock
            conda install conda-smithy
            conda-smithy rerender

    *  Get sha256 hash from pypi.org
    *  Update version number and hash in recipe
    *  Check dependencies
