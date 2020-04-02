FROM python:3.7

LABEL maintainer="Sebastian Salata R-T <sa.salatart@gmail.com>"

ENV PYTHON_ENV=production

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

COPY Pipfile Pipfile.lock ./

RUN apt-get update && pip3 install pipenv && pipenv install

COPY . .

CMD ["pipenv", "run", "start"]
