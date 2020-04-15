[mit]: https://img.shields.io/badge/License-MIT-blue.svg
[circleci]: https://circleci.com/gh/open-sir/open-sir.svg?style=shield
[black]: https://img.shields.io/badge/code%20style-black-000000.svg
[rtd]: https://readthedocs.org/projects/open-sir/badge/?version=latest

# open-sir

[![License: MIT][mit]](https://opensource.org/licenses/MIT)
[![open-sir][circleci]](https://circleci.com/gh/open-sir/open-sir)
[![Code style: black][black]](https://github.com/psf/black)
[![Documentation Status][rtd]](https://open-sir.readthedocs.io/en/latest)


Open-SIR is an Open Source Python project for modelling pandemics and
infectious diseases using Compartmental Models, such as the widely used
[Susceptible-Infected-Removed (SIR)
model](http://rocs.hu-berlin.de/corona/docs/forecast/model/#classic-sir-dynamics).
The current stage of the software is *Alpha*.

## Features
- Model the dynamics of infectious diseases
- Parameter fitting
- Calculation of confidence intervals
- CLI for interfacing with non Python environments (Bash, Node.JS, Matlab, etc).

So far, Open-SIR provides an implementation of the SIR model and the novel
[SIR-X model, developed by Maier and
Dirk](https://science.sciencemag.org/content/early/2020/04/07/science.abb4557.full)
from the [Robert Koch
Institut](http://rocs.hu-berlin.de/corona/docs/forecast/model/#sir-x-dynamics-outbreaks-with-temporally-increasing-interventions).

For the API reference and examples of usage, please check the
[online documentation](https://open-sir.readthedocs.io/).

### Getting started

Welcome to Open-Sir! Install Open-SIR using `pip`

```
git clone https://github.com/open-sir/open-sir.git
cd open-sir
pip install .
```

Please follow [Getting
Started](https://open-sir.readthedocs.io/en/latest/doc/getting-started.html) to
get your hands on Open-SIR.

## Authors

* **[José Álamos](https://github.com/jia200x)** -
  [RIOT-OS](https://github.com/RIOT-OS)
* **[Felipe Huerta](https://github.com/felipehuerta17)** - [PhD
  Student](https://www.imperial.ac.uk/people/f.huerta-perez17) at [Imperial
College London](https://github.com/ImperialCollegeLondon)
* **[Sebastián Salata](https://github.com/sasalatart)** - Software Engineer -
  Full Stack

See also the list of
[contributors](https://github.com/open-sir/open-sir/contributors) who
participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE)
file for details

## Acknowledgements

* [Robert Koch Institut](https://www.rki.de/EN/Home/homepage_node.html) for the
  clear explanation of SIR and SIR-X models.

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-4-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/jia200x"><img src="https://avatars3.githubusercontent.com/u/1260616?v=4" width="100px;" alt=""/><br /><sub><b>José Alamos</b></sub></a><br /><a href="https://github.com/open-sir/open-sir/commits?author=jia200x" title="Code">💻</a> <a href="https://github.com/open-sir/open-sir/commits?author=jia200x" title="Documentation">📖</a> <a href="#maintenance-jia200x" title="Maintenance">🚧</a></td>
    <td align="center"><a href="http://www.imperial.ac.uk/people/f.huerta-perez17"><img src="https://avatars3.githubusercontent.com/u/33637198?v=4" width="100px;" alt=""/><br /><sub><b>Felipe Huerta</b></sub></a><br /><a href="https://github.com/open-sir/open-sir/commits?author=felipehuerta17" title="Code">💻</a> <a href="https://github.com/open-sir/open-sir/commits?author=felipehuerta17" title="Documentation">📖</a> <a href="#tutorial-felipehuerta17" title="Tutorials">✅</a></td>
    <td align="center"><a href="https://github.com/sasalatart"><img src="https://avatars1.githubusercontent.com/u/5463900?v=4" width="100px;" alt=""/><br /><sub><b>Sebastián Salata</b></sub></a><br /><a href="https://github.com/open-sir/open-sir/commits?author=sasalatart" title="Code">💻</a> <a href="#maintenance-sasalatart" title="Maintenance">🚧</a><a href="#infra-sasalatart" title="Infra"> 🚇</a></td>
    <td align="center"><a href="https://github.com/leandrolanzieri"><img src="https://avatars1.githubusercontent.com/u/5381296?v=4" width="100px;" alt=""/><br /><sub><b>Leandro Lanzieri</b></sub></a><br /><a href="https://github.com/open-sir/open-sir/commits?author=leandrolanzieri" title="Code">💻</a> <a href="#maintenance-leandrolanzieri" title="Maintenance">🚧</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
