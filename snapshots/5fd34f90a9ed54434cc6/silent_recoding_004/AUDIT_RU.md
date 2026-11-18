# 004 — adversarial audit (15 пунктов)

1. **Finite-sample noise:** coverage и conditional power разделены. KL и
   exact CP имеют аналитические marginal гарантии; точная small-n PMF
   проверка независима от CI construction, но не заменяет общий proof.
2. **Density artifact:** размер2141 больше не входит в binomial confidence;
   полнота B/D/U и измеренная неизвестная масса остаются существенными.
3. **Event-size artifact:** фиксированный бюджет<=2 LC событий, first exit,
   stops и NO_EXIT сохранены. Нет подмены Hamming distance количеством событий.
4. **Hidden common input:** ancestry/clonal descendants нельзя считать IID.
   Общая ancestor выборка явно вложенная; union bound не требует её независимости.
5. **Mixture contexts:** наличие species/locus/repair/batch mixtures не закрыто;
   требуется допустимая law и удостоверенный event bias, не просто ещё больше n.
6. **Temporal nonstationarity:** fixed-n и одинаковая sampling law внутри start;
   drift/adaptive stopping не покрыты без отдельной гарантии.
7. **Binning/windowing:** B/BU не идентифицируются endpoint-only в общем случае.
   Пример003 остаётся препятствием; новые бинарные bins не стирают порядок.
8. **Leakage:** priority не использует expected usefulness; панели фиксируются
   до новых данных. Functional bins выбираются независимо от supply sample.
   Ретроспективная model-conditioned focal selection не является подтверждением
   на новых данных; discovery19 cases и планирование явно отделены.
9. **Multiplicity:** four supply стороны alpha_S/4, функциональная семья
   пересчитана1244/1344/1644. Все3 baseline сохранены. Post-data выбор CP/KL
   запрещён без новой совместной гарантии. Coverage!=power!=биологическая95%.
10. **Stronger baseline:** точная binomial CP существенно улучшает KL sufficient
    counts; это стандартная статистика. Полные TV bounds не объявлены нижней ценой.
11. **Latent explanation:** экономия является дизайном наблюдений, не causal
    доказательством механизма. Native competence и общий G не установлены.
12. **Carrier:** тот же focal protein/DNA, те же7 starting roles/контроли,
    три baseline. Истории, reads и protein objects считаются раздельно.
13. **Observable relevance:** непосредственно Phi2, loss/restore; rewrite
    выводится суммой. Более консервативные coarse границы допускают ложный HOLD,
    но не меняют желаемое утверждение при действительных предпосылках.
14. **Generalization:** числа относятся к одной фиксированной плановой law,
    будущей классификации и выбранному readout. Не универсальный laboratory n,
    не optimum по носителям и не новая теорема SFP.
15. **Falsification:** если class intervals не разрешают S93R или bias/readout
    не закрывает условия, соответствующий дизайн отвергается/остаётся HOLD.
    Неположительные coarse bounds не опровергают саму biological hypothesis.

Программные контроли: positive arithmetic example, null frequencies, false
ancestor nesting, unequal ancestor n, invalid/boolean/fractional/out-of-range
counts, отсутствующие события; оба метода fail closed. Контроли не являются
биологическими наблюдениями и не называют число assertions доказательством.

Исполнение: первый новый004 запуск остановился assert на ошибочном выводе
числа baseline из variant panel. No-op background находится только в starting
roles; family source исправлена до выпускного результата. Старого success
receipt при отказе не было. Следующий расчёт прошёл. После добавления CP count
interface выполнен один целевой новый replay004. Accepted R1/003 не исполнялись.

GO на статистический дизайн. Внешняя независимая приёмка004 ожидается.
HOLD biological/event feasibility. Нет отправки письма, покупки или вмешательства.
