---
layout: posts
title: "Embedding Gist into webpage Markdown"
date: 2019-01-25
excerpt: "How-to on imbedding Gist content into Minimal Mistakes markdown blog post."
---

Turns out it is pretty easy to imbed an existing Gist into a Minimal Mistakes blog post. Simply write a blog post as you normally would and where you want to embed the Gist, simply add the following markdown text. The long alpha-numeric string is unique to each Gist and can be retrieved from the Gist URL.

```
{% gist 2d7ca4eed82764833f8d93a6fea28f15 %}
```

Then the embedded Gist should appear like this:

{% gist 2d7ca4eed82764833f8d93a6fea28f15 %}
